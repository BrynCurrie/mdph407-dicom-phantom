#!/usr/bin/env python3
"""
MDPH407 — synthetic DICOM CT + RT Structure Set generator.

Run this on your own laptop or a lab PC that is *not* the Monaco workstation.
Copy the output folder onto the ER313 machine and import it (see the lab notes).

Presets
=======
    python make_phantom.py --preset cube
    python make_phantom.py --preset cylinder
    python make_phantom.py --config configs/cylinder_offset_ptv.json

What is written
===============
    <out>/CT_000.dcm ...     CT Image Storage, one file per slice
    <out>/RS.dcm             RT Structure Set Storage (External, PTV, OARs)
    <out>/manifest.json      Geometry and UIDs for the lab report

The tags this script sets are the ones Monaco actually uses on import:
PatientName / PatientID, FrameOfReferenceUID, Study/Series/SOP Instance UIDs,
ImagePositionPatient, ImageOrientationPatient, PixelSpacing, SliceThickness,
and the RTSTRUCT references that point every contour back at the matching CT
slice. Change a UID family independently and the overlay will break — that is
the lesson, not an inconvenience.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime
from typing import Any

import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import (
    CTImageStorage,
    ExplicitVRLittleEndian,
    RTStructureSetStorage,
    generate_uid,
)

PRESETS = {
    "cube": os.path.join("configs", "cube_trs398.json"),
    "cylinder": os.path.join("configs", "cylinder_offset_ptv.json"),
}

CT_SOP_CLASS = CTImageStorage
RS_SOP_CLASS = RTStructureSetStorage


def grid_axes(cfg: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    g = cfg["grid"]
    nx, ny, nz = int(g["nx"]), int(g["ny"]), int(g["nz"])
    ps, st = float(g["pixel_size_mm"]), float(g["slice_thickness_mm"])
    x = (np.arange(nx) - (nx - 1) / 2.0) * ps
    y = (np.arange(ny) - (ny - 1) / 2.0) * ps
    z = (np.arange(nz) - (nz - 1) / 2.0) * st
    return x, y, z


def image_position(x: np.ndarray, y: np.ndarray, z: np.ndarray, k: int) -> list[float]:
    return [float(x[0]), float(y[0]), float(z[k])]


def inside_cube(X, Y, Z, center, size) -> np.ndarray:
    cx, cy, cz = center
    sx, sy, sz = size
    return (
        (np.abs(X - cx) <= sx / 2.0)
        & (np.abs(Y - cy) <= sy / 2.0)
        & (np.abs(Z - cz) <= sz / 2.0)
    )


def inside_cylinder(X, Y, Z, center, radius, height, axis="z") -> np.ndarray:
    cx, cy, cz = center
    if axis == "z":
        r2 = (X - cx) ** 2 + (Y - cy) ** 2
        return (r2 <= radius**2) & (np.abs(Z - cz) <= height / 2.0)
    if axis == "y":
        r2 = (X - cx) ** 2 + (Z - cz) ** 2
        return (r2 <= radius**2) & (np.abs(Y - cy) <= height / 2.0)
    r2 = (Y - cy) ** 2 + (Z - cz) ** 2
    return (r2 <= radius**2) & (np.abs(X - cx) <= height / 2.0)


def inside_sphere(X, Y, Z, center, radius) -> np.ndarray:
    cx, cy, cz = center
    return (X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2 <= radius**2


def mask_of(shape: str, X, Y, Z, spec: dict) -> np.ndarray:
    c = spec.get("center_mm", [0.0, 0.0, 0.0])
    if shape == "cube":
        return inside_cube(X, Y, Z, c, spec["size_mm"])
    if shape == "cylinder":
        return inside_cylinder(
            X, Y, Z, c, spec["radius_mm"], spec["height_mm"], spec.get("axis", "z")
        )
    if shape == "sphere":
        return inside_sphere(X, Y, Z, c, spec["radius_mm"])
    raise ValueError(f"Unknown shape: {shape}")


def circle_points(cx, cy, z, radius, n=72) -> list[float]:
    pts: list[float] = []
    for i in range(n):
        a = 2.0 * math.pi * i / n
        pts.extend([cx + radius * math.cos(a), cy + radius * math.sin(a), z])
    return pts


def square_points(cx, cy, z, sx, sy) -> list[float]:
    hx, hy = sx / 2.0, sy / 2.0
    corners = [
        (cx - hx, cy - hy),
        (cx + hx, cy - hy),
        (cx + hx, cy + hy),
        (cx - hx, cy + hy),
    ]
    pts: list[float] = []
    for x, y in corners:
        pts.extend([x, y, z])
    return pts


def contours_for(spec: dict, z_axis: np.ndarray) -> list[tuple[int, list[float]]]:
    shape = spec["shape"]
    cx, cy, cz = spec.get("center_mm", [0.0, 0.0, 0.0])
    out: list[tuple[int, list[float]]] = []
    for k, z in enumerate(z_axis):
        if shape == "cube":
            sx, sy, sz = spec["size_mm"]
            if abs(z - cz) > sz / 2.0:
                continue
            out.append((k, square_points(cx, cy, float(z), sx, sy)))
        elif shape == "cylinder" and spec.get("axis", "z") == "z":
            if abs(z - cz) > spec["height_mm"] / 2.0:
                continue
            out.append((k, circle_points(cx, cy, float(z), spec["radius_mm"])))
        elif shape == "sphere":
            r = spec["radius_mm"]
            half = r * r - (z - cz) ** 2
            if half <= 0:
                continue
            out.append((k, circle_points(cx, cy, float(z), math.sqrt(half))))
        else:
            raise ValueError(f"Contouring not implemented for {shape}")
    return out


def build_volume(cfg: dict, x, y, z) -> np.ndarray:
    X, Y, Z = np.meshgrid(x, y, z, indexing="xy")
    vol = np.full(X.shape, int(cfg.get("background_hu", -1000)), dtype=np.int16)
    body = cfg["body"]
    vol[mask_of(body["shape"], X, Y, Z, body)] = int(body.get("hu", 0))
    for spec in cfg.get("structures", []):
        vol[mask_of(spec["shape"], X, Y, Z, spec)] = int(spec.get("hu", 0))
    return vol


def _file_meta(sop_class: str, sop_instance: str) -> FileMetaDataset:
    meta = FileMetaDataset()
    meta.FileMetaInformationGroupLength = 0
    meta.FileMetaInformationVersion = b"\x00\x01"
    meta.MediaStorageSOPClassUID = sop_class
    meta.MediaStorageSOPInstanceUID = sop_instance
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    meta.ImplementationClassUID = generate_uid()
    meta.ImplementationVersionName = "MDPH407_1.0"
    return meta


def write_ct_slices(cfg, vol, x, y, z, out_dir, study_uid, series_uid, for_uid, now):
    g = cfg["grid"]
    ps, st = float(g["pixel_size_mm"]), float(g["slice_thickness_mm"])
    ny, nx, nz = vol.shape
    sop_uids = []
    for k in range(nz):
        sop = generate_uid()
        sop_uids.append(sop)
        meta = _file_meta(CT_SOP_CLASS, sop)
        ds = FileDataset(None, {}, file_meta=meta, preamble=b"\x00" * 128)
        ds.is_little_endian = True
        ds.is_implicit_VR = False
        ds.SOPClassUID = CT_SOP_CLASS
        ds.SOPInstanceUID = sop
        ds.StudyInstanceUID = study_uid
        ds.SeriesInstanceUID = series_uid
        ds.FrameOfReferenceUID = for_uid
        ds.SpecificCharacterSet = "ISO_IR 100"
        ds.Modality = "CT"
        ds.Manufacturer = "University of Canterbury"
        ds.InstitutionName = "MDPH407 educational"
        ds.StudyDescription = "MDPH407 synthetic phantom"
        ds.SeriesDescription = cfg.get("description", "synthetic CT")
        ds.StudyDate = now.strftime("%Y%m%d")
        ds.SeriesDate = ds.StudyDate
        ds.ContentDate = ds.StudyDate
        ds.StudyTime = now.strftime("%H%M%S")
        ds.SeriesTime = ds.StudyTime
        ds.ContentTime = ds.StudyTime
        ds.AccessionNumber = ""
        ds.ReferringPhysicianName = ""
        ds.StudyID = "MDPH407"
        ds.SeriesNumber = 1
        ds.InstanceNumber = k + 1
        ds.AcquisitionNumber = 1
        ds.PatientName = cfg.get("patient_name", "MDPH407^PHANTOM")
        ds.PatientID = cfg.get("patient_id", "MDPH407PHANTOM")
        ds.PatientBirthDate = ""
        ds.PatientSex = "O"
        ds.PatientPosition = "HFS"
        ds.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.Rows = ny
        ds.Columns = nx
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1
        ds.RescaleIntercept = 0.0
        ds.RescaleSlope = 1.0
        ds.WindowCenter = 40
        ds.WindowWidth = 400
        ds.PixelSpacing = [ps, ps]
        ds.SliceThickness = st
        ds.SpacingBetweenSlices = st
        ds.SliceLocation = float(z[k])
        ds.ImagePositionPatient = image_position(x, y, z, k)
        ds.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        ds.PositionReferenceIndicator = ""
        slc = np.ascontiguousarray(vol[:, :, k], dtype=np.int16)
        ds.PixelData = slc.tobytes()
        ds.save_as(os.path.join(out_dir, f"CT_{k:03d}.dcm"), write_like_original=False)
    return sop_uids


def _roi_observation(number, name, rt_type):
    obs = Dataset()
    obs.ObservationNumber = number
    obs.ReferencedROINumber = number
    obs.ROIObservationLabel = name
    obs.RTROIInterpretedType = rt_type
    obs.ROIInterpreter = ""
    return obs


def write_rtstruct(cfg, x, y, z, out_dir, study_uid, ct_series_uid, rs_series_uid, for_uid, sop_uids, now):
    sop = generate_uid()
    meta = _file_meta(RS_SOP_CLASS, sop)
    ds = FileDataset(None, {}, file_meta=meta, preamble=b"\x00" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.SOPClassUID = RS_SOP_CLASS
    ds.SOPInstanceUID = sop
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = rs_series_uid
    ds.FrameOfReferenceUID = for_uid
    ds.Modality = "RTSTRUCT"
    ds.Manufacturer = "University of Canterbury"
    ds.InstitutionName = "MDPH407 educational"
    ds.StructureSetLabel = cfg.get("studyset_hint", "PHANTOM")[:16]
    ds.StructureSetName = cfg.get("preset", "phantom")
    ds.StructureSetDate = now.strftime("%Y%m%d")
    ds.StructureSetTime = now.strftime("%H%M%S")
    ds.ApprovalStatus = "UNAPPROVED"
    ds.PatientName = cfg.get("patient_name", "MDPH407^PHANTOM")
    ds.PatientID = cfg.get("patient_id", "MDPH407PHANTOM")
    ds.PatientBirthDate = ""
    ds.PatientSex = "O"
    images = []
    for u in sop_uids:
        ci = Dataset()
        ci.ReferencedSOPClassUID = CT_SOP_CLASS
        ci.ReferencedSOPInstanceUID = u
        images.append(ci)
    ref_series = Dataset()
    ref_series.SeriesInstanceUID = ct_series_uid
    ref_series.ContourImageSequence = Sequence(images)
    ref_study = Dataset()
    ref_study.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.1"
    ref_study.ReferencedSOPInstanceUID = study_uid
    ref_study.RTReferencedSeriesSequence = Sequence([ref_series])
    ref_for = Dataset()
    ref_for.FrameOfReferenceUID = for_uid
    ref_for.RTReferencedStudySequence = Sequence([ref_study])
    ds.ReferencedFrameOfReferenceSequence = Sequence([ref_for])
    rois = [cfg["body"]] + list(cfg.get("structures", []))
    ss_roi, roi_contour, observations = [], [], []
    for i, spec in enumerate(rois, start=1):
        roi = Dataset()
        roi.ROINumber = i
        roi.ReferencedFrameOfReferenceUID = for_uid
        roi.ROIName = spec["name"]
        roi.ROIGenerationAlgorithm = "AUTOMATIC"
        roi.ROIDescription = spec.get("shape", "")
        ss_roi.append(roi)
        rc = Dataset()
        rc.ReferencedROINumber = i
        rc.ROIDisplayColor = [int(c) for c in spec.get("color", [255, 255, 0])]
        contour_seq = []
        for k, xyz in contours_for(spec, z):
            c = Dataset()
            c.ContourGeometricType = "CLOSED_PLANAR"
            c.NumberOfContourPoints = len(xyz) // 3
            c.ContourData = [float(v) for v in xyz]
            img = Dataset()
            img.ReferencedSOPClassUID = CT_SOP_CLASS
            img.ReferencedSOPInstanceUID = sop_uids[k]
            c.ContourImageSequence = Sequence([img])
            contour_seq.append(c)
        rc.ContourSequence = Sequence(contour_seq)
        roi_contour.append(rc)
        observations.append(_roi_observation(i, spec["name"], spec.get("rt_type", "ORGAN")))
    ds.StructureSetROISequence = Sequence(ss_roi)
    ds.ROIContourSequence = Sequence(roi_contour)
    ds.RTROIObservationsSequence = Sequence(observations)
    ds.save_as(os.path.join(out_dir, "RS.dcm"), write_like_original=False)
    return sop


def write_manifest(cfg, out_dir, uids, now):
    doc = {
        "generated": now.isoformat(timespec="seconds"),
        "script": "make_phantom.py",
        "preset": cfg.get("preset"),
        "description": cfg.get("description"),
        "patient_name": cfg.get("patient_name"),
        "patient_id": cfg.get("patient_id"),
        "grid": cfg["grid"],
        "body": cfg["body"],
        "structures": cfg.get("structures", []),
        "uids": uids,
        "notes": [
            "PatientName and PatientID must match between CT and RTSTRUCT.",
            "FrameOfReferenceUID must match or Monaco will not overlay contours.",
            "StudyInstanceUID is shared; SeriesInstanceUID differs for CT vs RS.",
            "Do not reuse these UIDs for a second phantom — mint new ones.",
        ],
    }
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2)


def load_config(args):
    if args.config:
        path = args.config
    elif args.preset:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(here, PRESETS[args.preset])
    else:
        raise SystemExit("Pass --preset cube|cylinder or --config path.json")
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    if args.patient_id:
        cfg["patient_id"] = args.patient_id
    if args.patient_name:
        cfg["patient_name"] = args.patient_name
    return cfg


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--preset", choices=sorted(PRESETS), help="Built-in lab phantom")
    p.add_argument("--config", help="JSON config (shape, target, OARs)")
    p.add_argument("--out", default=None, help="Output folder (default: output/<preset>)")
    p.add_argument("--patient-id", help="Override PatientID (use your student ID)")
    p.add_argument("--patient-name", help="Override PatientName, DICOM style GIVEN^FAMILY")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args)
    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out or os.path.join(here, "output", cfg.get("preset", "phantom"))
    os.makedirs(out, exist_ok=True)
    x, y, z = grid_axes(cfg)
    vol = build_volume(cfg, x, y, z)
    now = datetime.now()
    study_uid = generate_uid()
    ct_series_uid = generate_uid()
    rs_series_uid = generate_uid()
    for_uid = generate_uid()
    sop_uids = write_ct_slices(cfg, vol, x, y, z, out, study_uid, ct_series_uid, for_uid, now)
    rs_sop = write_rtstruct(cfg, x, y, z, out, study_uid, ct_series_uid, rs_series_uid, for_uid, sop_uids, now)
    write_manifest(cfg, out, {
        "StudyInstanceUID": study_uid,
        "CT_SeriesInstanceUID": ct_series_uid,
        "RS_SeriesInstanceUID": rs_series_uid,
        "FrameOfReferenceUID": for_uid,
        "RS_SOPInstanceUID": rs_sop,
    }, now)
    print(f"Wrote {len(sop_uids)} CT slices + RS.dcm → {out}")
    print(f"PatientID={cfg.get('patient_id')}  FoR={for_uid}")


if __name__ == "__main__":
    main()
