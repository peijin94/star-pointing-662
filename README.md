# EOVSA ASI662 Star Pointing

This project uses a ZWO ASI662MC camera and a short-focal-length lens to measure pointing errors for an Expanded Owens Valley Solar Array (EOVSA) antenna. It captures scheduled star fields, solves each image with ASTAP, compares the commanded and measured sky positions, converts those differences to altitude/azimuth offsets, and fits a seven-term pointing model.

## Status

The original code was copied from the Raspberry Pi development directory at `/home/starp/star_capture` on 2026-08-20. Hardware-free regression tests run in the Pi's Python environment, but the corrected acquisition and pointing workflow has not yet been exercised with the camera, ASTAP, or an antenna. See [Validation still required](#validation-still-required-before-antenna-use).

The code measures pointing; it does **not** command or slew the EOVSA antenna. Antenna scheduling/control must provide the target table and move the antenna independently.

## Workflow

1. Read a scheduled star-target table.
2. Wait for the antenna settling window associated with each target.
3. Capture a 16-bit, 1920 x 1080 image from the ASI662MC.
4. Write the image and approximate target coordinates to FITS.
5. Run ASTAP and read the solved center coordinates from its WCS output.
6. Record nominal and solved right ascension/declination.
7. Convert the coordinate differences to altitude/azimuth pointing offsets.
8. Fit an IA/IE/CA/NPAE/AN/AW/TF pointing model.

## Repository contents

- `starpoint_662.py` — target-table parsing, image capture, FITS creation, ASTAP solving, coordinate conversion, and pointing-model fitting.
- `focus_662.py` — continuously captures and displays images for camera focusing.
- `test_focus_662.py` — regression checks for focus-script CLI argument binding and gain conversion.
- `test_starpoint_662.py` — hardware-free regression coverage for parsing, camera lifecycle, FITS metadata, ASTAP handling, output paths, coordinate wrapping, and pointing-model fitting.

## Hardware and external software

- Raspberry Pi with the EOVSA camera interface
- ZWO ASI662MC camera
- Lens currently represented as 25 mm focal length in the FITS headers
- ZWO ASI SDK, expected at `/usr/local/lib/libASICamera2.so`
- ASTAP executable available as `astap` on `PATH`
- An ASTAP star database appropriate for the camera field of view
- X11 display forwarding if plots or the focus window will be shown remotely

Python dependencies used by the source are:

```text
zwoasi
numpy
astropy
matplotlib
scipy
```

Use the camera host's existing virtual environment where possible. The current Pi setup has used:

```bash
source /home/dgary/star_capture/starenv/bin/activate
```

## Focusing

From the repository directory on the Pi:

```bash
python focus_662.py --exp 0.01 --gain 252 --scale log
```

The focus program captures continuously. Stop it with Ctrl-C.

## Star-pointing functions

`starpoint_662.py` currently exposes Python functions rather than a command-line entry point. A typical interactive workflow is:

```python
from starpoint_662 import imaging, result2dazel, mountcal

imaging("startable-YYYY-MM-DD.txt", exp_s=3, test=True)
```

Use `test=True` first to exercise schedule parsing and timing without opening the camera or running ASTAP. A hardware acquisition uses:

```python
imaging("startable-YYYY-MM-DD.txt", exp_s=3, test=False)
```

After acquisition from the repository root, the intended post-processing calls are:

```python
result_file = "YYYY-MM-DD/YYYY-MM-DD_solve_results.txt"
result2dazel(result_file)
coefficients = mountcal(result_file)
```

`imaging()` writes through explicit paths and does not change the process working directory.

## Target-table contract

`startable2dict()` expects a fixed-width text format:

- The first two lines are ignored.
- The first four characters of each later line form the target name.
- Characters 17 through 40 contain whitespace-separated `HH:MM:SS` right ascension and signed `DD:MM:SS` declination.
- The final two whitespace-separated fields form an Astropy-compatible UTC date and time.

The parser validates the coordinate ranges and reports malformed rows with their filename and line number. Preserve the fixed columns and verify the parsed targets before acquisition.

## Outputs

During acquisition, the code creates a directory named for the first target date without changing the process working directory. It may then create:

- FITS images containing camera and approximate pointing metadata
- ASTAP `.wcs` solution files
- `<date>_solve_results.txt` containing nominal and solved sky coordinates
- `<name>_dazel.txt` containing calculated altitude/azimuth offsets
- A two-panel Matplotlib calibration plot

The observing site, atmosphere, camera, and lens assumptions are currently constants in `starpoint_662.py`. Confirm `LAT`, `LON`, `TEMP_C`, `PRESS_MB`, focal length, pixel size, image dimensions, gain, exposure, and ASTAP field of view before an observing run.

## Corrected review findings

The correction pass added guaranteed camera cleanup, zero-coordinate FITS support, accurate gain metadata, safe logarithmic display, ASTAP timeout/exit/stale-file handling, contextual table validation, stable output paths, first-row preservation, wrapped azimuth differences, true arcminute reporting, minimum-data/convergence/rank checks, and observed-versus-model plots.

## Validation still required before antenna use

1. **Pointing conventions:** verify the seven model-term signs, axis conventions, and coefficient definitions against the EOVSA control system before applying coefficients.
2. **Refraction:** confirm whether refraction belongs in the nominal altitude used by this model and whether pressure, temperature, and site elevation should be supplied through Astropy's `AltAz` frame instead.
3. **Site and optics:** confirm latitude, longitude, site height, focal length, pixel size, image dimensions, gain, exposure, ASTAP field of view, and search radius against the installed antenna/camera assembly.
4. **Schedule semantics:** confirm that each table timestamp denotes the start of antenna settling and that acquisition should occur only from 30 to 60 seconds afterward.
5. **External integration:** run controlled checks with the actual ASI662MC and ASTAP database. The automated tests mock the camera and solver and do not slew the antenna.
6. **Fit acceptance:** convergence and rank are checked, but operational acceptance thresholds for residual RMS, coefficient uncertainty, sky coverage, and outlier rejection still need domain decisions.

Do not apply fitted coefficients to an operational antenna until these items are validated.

## Current test coverage

Run all tests in the Pi environment without opening the camera:

```bash
MPLBACKEND=Agg python -m unittest -v test_focus_662.py test_starpoint_662.py
```

The suite covers:

- exposure and gain binding to the correct function parameters
- parsing camera gain as an integer required by the ZWO SDK
- camera closure on acquisition failures
- schedule and solve-result validation
- zero-valued FITS coordinates and requested gain metadata
- stale and successful ASTAP WCS handling
- stable process/output paths
- first-observation preservation and azimuth wraparound
- radians-to-arcminutes conversion
- minimum fit data and recovery of seven synthetic pointing coefficients

The suite does not operate the camera, invoke a real ASTAP process, or command an antenna.
