# MDPH407 DICOM phantom

Teaching generator for the University of Canterbury **MDPH407** treatment-planning lab (Elekta Monaco 6.2.3.0, educational licence).

It writes a synthetic **CT Image** series and an **RT Structure Set** that share one Frame of Reference. You choose a body shape, a target, and OARs in a JSON config. The lab notes tell you how to import the folder into Monaco.

Do not run this *on* the Monaco PC if you can avoid it. Generate the files on your laptop, copy the folder to ER313, import.

Repository: https://github.com/BrynCurrie/mdph407-dicom-phantom

## Why this exists

Monaco will not calculate dose without an External contour sitting on a CT that has a CT-to-ED curve. A clinical CT already has both. This script is the stand-in: it mints legal DICOM so you can spend lab time on beams rather than on drawing a 40 cm cube by hand — and so you can *see* the tags that have to agree.

| Object | File | SOP class |
| --- | --- | --- |
| CT slices | `CT_000.dcm` … | CT Image Storage |
| Structures | `RS.dcm` | RT Structure Set Storage |
| Geometry + UIDs | `manifest.json` | (not DICOM; for the report) |

The four DICOM RT objects of a finished plan are CT, RTSTRUCT, RTPLAN, RTDOSE. This repo produces the first two. Monaco produces the last two after you plan.

## Install

Python 3.10+ on your own machine.

```bash
git clone https://github.com/BrynCurrie/mdph407-dicom-phantom.git
cd mdph407-dicom-phantom
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Lab presets

```bash
python make_phantom.py --preset cube     --patient-id MDPH407123 --out output/cube
python make_phantom.py --preset cylinder --patient-id MDPH407123 --out output/cyl
```

Use your own student ID as `--patient-id` so the nine ER313 databases do not collide.

| Preset | Body | Target | OARs | Use in the lab |
| --- | --- | --- | --- | --- |
| `cube` | 40 cm water cube, HU 0 | none | none | TRS-398 10×10 cm field, PDD, SSD+10 cm, parallel opposed pair |
| `cylinder` | 40 cm ∅ × 40 cm water cylinder | 10 cm radius sphere, +5 cm in x | `OAR_Cord` (10 mm ∅ cylinder) and `OAR_Lung` (40 mm sphere, HU −700) | Cover the PTV, spare the OARs, compare with the class |

Both presets: 2 mm pixels, 5 mm slices, 200×200×80, HFS, air HU = −1000 outside the body. Monaco needs at least five slices; these have eighty.

## Your own shape / target / OARs

Copy a file in `configs/` and point `--config` at it. Supported shapes: `cube`, `cylinder`, `sphere`.

`rt_type` becomes `RTROIInterpretedType` on the structure set. Use `EXTERNAL` for the body, `PTV` for the target, `ORGAN` (or `AVOIDANCE`) for OARs.

Coordinates are DICOM patient millimetres, HFS: **+x right, +y posterior, +z superior**, origin at the centre of the phantom.

## Tags that must agree

| Tag | Why it matters |
| --- | --- |
| PatientName, PatientID | Must be identical on every CT slice and on `RS.dcm` or Monaco treats them as different patients |
| FrameOfReferenceUID | Shared by CT and RTSTRUCT. Different FoR → contours will not overlay |
| StudyInstanceUID | Shared. SeriesInstanceUID is *different* for the CT series and the structure set |
| SOPInstanceUID | Unique per file. The structure set references each slice SOP it contours |
| ImageOrientationPatient `[1,0,0,0,1,0]` + PatientPosition `HFS` | Stops left/right flipping in BEV |
| PixelSpacing, SliceThickness, ImagePositionPatient | Build the patient coordinate system the contours are drawn in |

The script mints new UIDs on every run (`pydicom.uid.generate_uid`). Do not copy UIDs out of a hospital study.

## Import into Monaco 6.2.3.0

1. Copy the output folder onto the ER313 PC you are sitting at (MAC-locked licence — files do not follow you to the next desk).
2. Open Patient → Import New Data.
3. Select the folder. Add the **CT series and the RT Structure Set** to the Transfer list.
4. Clinic / Patient ID / studyset name (≤14 characters; `studyset_hint` in the JSON is a suggestion).
5. Assign the CT-to-ED file the instructor names.
6. Import, open, confirm External / PTV / OARs overlay the slices.

If import fails, Show Log. The usual faults are a missing FoR, a Patient ID clash with something already in that PC’s database, or a series with fewer than five slices.

After import you still Force ED = 1.000 on External (and on any structure you want to treat as water). The lung OAR on the cylinder preset is HU −700 so you can see whether the algorithm cares.

## What this is not

- Not a clinical phantom.
- Not a substitute for commissioning a CT-to-ED curve.
- Not a licence to skip contour review. Look at the ends of the stack.

Product documentation for Monaco remains © Elekta and lives on Learn. This repo is course material.
