# EOVSA ASI662 Star Pointing

This project uses a ZWO ASI662MC camera and a short-focal-length lens to measure pointing errors for an Expanded Owens Valley Solar Array (EOVSA) antenna. It captures scheduled star fields, solves each image with ASTAP, compares the commanded and measured sky positions, converts those differences to altitude/azimuth offsets, and fits a seven-term pointing model.

## Status

The current code was copied from the Raspberry Pi development directory at `/home/starp/star_capture` on 2026-08-20. The source has received a static review, but `starpoint_662.py` has not yet been executed from this repository or validated on an antenna. See [Known issues](#known-issues-before-hardware-use) before running it.

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

After acquisition, the intended post-processing calls in the **same Python process** are:

```python
result_file = "YYYY-MM-DD_solve_results.txt"
result2dazel(result_file)
coefficients = mountcal(result_file)
```

This works because `imaging()` has already changed into the date directory. From a fresh process started at the repository root, use `YYYY-MM-DD/YYYY-MM-DD_solve_results.txt` instead. Review the working-directory issue below before using either form.

## Target-table contract

`startable2dict()` expects a fixed-width text format:

- The first two lines are ignored.
- The first four characters of each later line form the target name.
- Characters 17 through 40 contain whitespace-separated `HH:MM:SS` right ascension and signed `DD:MM:SS` declination.
- The final two whitespace-separated fields form an Astropy-compatible UTC date and time.

The parser does not currently validate this structure. Preserve the fixed columns and verify the parsed targets before acquisition.

## Outputs

During acquisition, the code creates a directory named for the first target date and changes into it. It may then create:

- FITS images containing camera and approximate pointing metadata
- ASTAP `.wcs` solution files
- `<date>_solve_results.txt` containing nominal and solved sky coordinates
- `<name>_dazel.txt` containing calculated altitude/azimuth offsets
- A two-panel Matplotlib calibration plot

The observing site, atmosphere, camera, and lens assumptions are currently constants in `starpoint_662.py`. Confirm `LAT`, `LON`, `TEMP_C`, `PRESS_MB`, focal length, pixel size, image dimensions, gain, exposure, and ASTAP field of view before an observing run.

## Known issues before hardware use

The following were found by static review and have not yet been corrected:

1. **Camera cleanup:** `snap()` does not close the camera in a `finally` block. A capture, reshape, FITS-write, or plotting exception can leave the camera open.
2. **Zero-valued coordinates:** checks such as `if ra` and `if dec` treat valid zero-degree coordinates as missing.
3. **Incorrect FITS gain metadata:** the `GAIN` header is always written as 252 even when another gain is passed to `snap()`.
4. **Stale plate solutions:** `plate_solve()` does not check the ASTAP exit status, remove an old WCS file before solving, or impose a timeout. A stale `.wcs` file could be mistaken for a new successful solution.
5. **Working-directory side effect:** `imaging()` changes into its date output directory and never restores the caller's original directory. Repeated calls in one process can use unexpected nested paths.
6. **Empty or malformed schedules:** an empty target list fails at `targets[0]`, while malformed fixed-width rows fail without a contextual error message.
7. **First calibration row discarded:** acquisition writes no results header, but `get_alt_az_data()` always skips the first line as if it were a header.
8. **Azimuth wraparound:** direct subtraction of azimuths near 0/360 degrees can create an error close to a full revolution instead of the small wrapped difference.
9. **Incorrect displayed units:** `mountcal()` labels coefficients as arcminutes but converts radians only to degrees; arcminutes require an additional factor of 60.
10. **Fit validation:** the seven-term least-squares fit does not enforce enough observations, report conditioning/uncertainty, or verify convergence and residual quality.
11. **Refraction/model validation:** refraction is applied asymmetrically to the model altitude, and the pointing-equation signs and conventions have not yet been validated against the EOVSA control system.
12. **Log display:** `np.log10()` receives raw image values directly; zero-valued pixels produce negative infinity and warnings.

These issues should be resolved and covered by tests before applying fitted coefficients to an operational antenna.

## Current test coverage

`test_focus_662.py` covers two previously fixed focus-script failures:

- exposure and gain binding to the correct function parameters
- parsing camera gain as an integer required by the ZWO SDK

There are not yet tests for target parsing, timing windows, FITS headers, ASTAP failure behavior, coordinate wraparound, or pointing-model recovery from synthetic data.
