import shutil
import subprocess
import sys
from pathlib import Path
from time import sleep

import astropy.units as u
import matplotlib.pylab as plt
import numpy as np
import zwoasi as asi
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time, TimeDelta
from scipy.optimize import leastsq


ASI_LIBRARY_PATH = "/usr/local/lib/libASICamera2.so"
IMAGE_SHAPE = (1080, 1920)
ASTAP_TIMEOUT_S = 120

# Observing-site and atmospheric configuration.
LAT = 37.5177
LON = -118.3655
TEMP_C = 10.0
PRESS_MB = 870


def print_dual(left_text, right_text):
    """Print left- and right-aligned status text on one terminal line."""
    columns, _ = shutil.get_terminal_size()
    padding = max(columns - len(left_text) - len(right_text), 1)
    sys.stdout.write(f"\r{left_text}{' ' * padding}{right_text}")
    sys.stdout.flush()


def _parse_sexagesimal(value, label):
    parts = value.split(":")
    if len(parts) != 3:
        raise ValueError(f"{label} must have three colon-separated fields")
    return [float(part) for part in parts]


def startable2dict(filename):
    """Read a fixed-width star table and return its scheduled targets."""
    path = Path(filename)
    targets = []

    with path.open("r", encoding="utf-8") as stream:
        lines = stream.readlines()

    for line_number, line in enumerate(lines[2:], start=3):
        if not line.strip():
            continue

        try:
            fields = line.split()
            if len(fields) < 2:
                raise ValueError("missing UTC date/time")

            name = line[:4].strip()
            if not name:
                raise ValueError("missing target name")

            coordinates = line[17:41].split()
            if len(coordinates) != 2:
                raise ValueError("missing fixed-width RA/Dec fields")

            rastr, decstr = coordinates
            ra_h, ra_m, ra_s = _parse_sexagesimal(rastr, "RA")
            dec_d_signed, dec_m, dec_s = _parse_sexagesimal(decstr, "Dec")
            dec_negative = decstr.startswith("-")
            dec_d = abs(dec_d_signed)

            if not (0 <= ra_h < 24 and 0 <= ra_m < 60 and 0 <= ra_s < 60):
                raise ValueError("RA is outside its valid range")
            if not (0 <= dec_d <= 90 and 0 <= dec_m < 60 and 0 <= dec_s < 60):
                raise ValueError("Dec is outside its valid range")
            if dec_d == 90 and (dec_m != 0 or dec_s != 0):
                raise ValueError("Dec exceeds 90 degrees")

            target_time = Time(f"{fields[-2]} {fields[-1]}")
        except (IndexError, TypeError, ValueError) as error:
            raise ValueError(
                f"{path}:{line_number}: malformed target row: {error}"
            ) from error

        targets.append(
            {
                "name": name,
                "time": target_time,
                "ra_h": int(ra_h),
                "ra_m": int(ra_m),
                "ra_s": ra_s,
                "dec_neg": dec_negative,
                "dec_d": int(dec_d),
                "dec_m": int(dec_m),
                "dec_s": dec_s,
            }
        )

    if not targets:
        raise ValueError(f"{path}: no targets found")

    return targets


def target_coordinates_deg(target):
    """Convert a parsed target's sexagesimal coordinates to degrees."""
    ra = target["ra_h"] * 15 + target["ra_m"] / 4.0 + target["ra_s"] / 240.0
    dec = target["dec_d"] + target["dec_m"] / 60.0 + target["dec_s"] / 3600.0
    if target["dec_neg"]:
        dec = -dec
    return ra, dec


def init_camera():
    """Initialize and return the first connected ZWO ASI camera."""
    asi.init(ASI_LIBRARY_PATH)
    if asi.get_num_cameras() == 0:
        raise ValueError("No cameras found")
    return asi.Camera(0)


def _capture_filename(srcname, ra, dec, exp_s):
    source = Path(srcname)
    exposure = f"{exp_s:g}"
    if ra is not None and dec is not None:
        suffix = f"_{int(ra)}_{int(dec)}_{exposure}.fits"
    else:
        suffix = f"_unk_unk_{exposure}.fits"
    return source.with_name(source.name + suffix)


def snap(srcname="test", ra=None, dec=None, gain=252, exp_s=3, show=True):
    """Capture one ASI662MC image and save it as FITS."""
    camera = init_camera()
    try:
        camera.set_control_value(asi.ASI_GAIN, gain)
        camera.set_control_value(asi.ASI_EXPOSURE, int(exp_s * 1_000_000))
        camera.set_image_type(asi.ASI_IMG_RAW16)

        print("Exposing", exp_s, "seconds on", srcname)
        buffer = camera.capture()
        image_data = np.frombuffer(buffer, dtype=np.uint16).reshape(IMAGE_SHAPE)

        hdu = fits.PrimaryHDU(image_data)
        header = hdu.header
        header["FOCALLEN"] = 25.0
        if ra is not None:
            header["RA"] = ra
        if dec is not None:
            header["DEC"] = dec
        header["INSTRUME"] = "ASI662MC"
        header["GAIN"] = gain
        header["EXPTIME"] = exp_s
        header["XPIXSZ"] = 2.9
        header["YPIXSZ"] = 2.9

        filename = _capture_filename(srcname, ra, dec, exp_s)
        print("Writing file", filename)
        hdu.writeto(filename, overwrite=True)

        if show:
            plt.ion()
            plt.imshow(
                np.log10(np.maximum(image_data, 1)),
                origin="lower",
                cmap="gray",
            )
            plt.pause(1)

        return str(filename)
    finally:
        camera.close()


def _read_wcs_center(wcs_file):
    solution = {}
    with Path(wcs_file).open("r", encoding="utf-8") as stream:
        for line in stream:
            key, separator, remainder = line.partition("=")
            if separator:
                solution[key.strip()] = remainder.partition("/")[0].strip()

    try:
        return float(solution["CRVAL1"]), float(solution["CRVAL2"])
    except (KeyError, TypeError, ValueError):
        return None


def plate_solve(filename, timeout_s=ASTAP_TIMEOUT_S):
    """Run ASTAP and return the solved center, or ``None`` on solve failure."""
    fits_file = Path(filename)
    wcs_file = fits_file.with_suffix(".wcs")
    wcs_file.unlink(missing_ok=True)
    command = ["astap", "-f", str(fits_file), "-fov", "7.2", "-r", "3"]

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None

    if result.returncode != 0 or not wcs_file.is_file():
        return None

    return _read_wcs_center(wcs_file)


def imaging(startable, exp_s=3, test=False):
    """Run scheduled capture and plate solving for every target in a star table."""
    targets = startable2dict(startable)
    settle_time = TimeDelta(30 * u.s)
    date_string = targets[0]["time"].iso.split()[0]
    output_folder = Path(date_string)
    output_folder.mkdir(exist_ok=True)
    results_file = output_folder / f"{date_string}_solve_results.txt"

    with results_file.open("w", encoding="utf-8") as stream:
        for target in targets:
            while Time.now() < target["time"] + settle_time:
                sleep(1)
                print_dual(
                    f"Waiting for {target['name']} time {target['time'].iso[:19]}",
                    Time.now().iso[:19],
                )

            if Time.now() >= target["time"] + 2 * settle_time:
                print("Skipped", target["name"], target["time"].iso[:19])
                continue

            ra, dec = target_coordinates_deg(target)
            if test:
                print(f"{target['name']}: {ra:8.4f} {dec:8.4f} test")
            else:
                time_string = target["time"].iso.split()[1].replace(":", "")
                source_name = output_folder / f"{target['name']}_{time_string[:4]}"
                fits_file = snap(
                    srcname=str(source_name),
                    ra=ra,
                    dec=dec,
                    gain=252,
                    exp_s=exp_s,
                    show=False,
                )
                result = plate_solve(fits_file)
                if result is not None:
                    actual_ra, actual_dec = result
                    success_line = (
                        f"{target['name']}: {ra:8.4f} {dec:8.4f} actual: "
                        f"{actual_ra:8.4f} {actual_dec:8.4f} "
                        f"{target['time'].iso[:19]}"
                    )
                    print(success_line)
                    stream.write(success_line + "\n")
                    stream.flush()
                else:
                    print(target["name"], "solve failed!")

            while Time.now() < target["time"] + 2 * settle_time:
                sleep(1)
                print_dual(
                    f"Completed {target['name']} time {target['time'].iso[:19]}",
                    Time.now().iso[:19],
                )
            print()

    print("All Done!")
    return str(results_file)


def get_refraction_deg(alt_deg):
    """Return Bennett atmospheric refraction in degrees."""
    if alt_deg <= 0:
        return 0
    refraction_arcmin = 1.0 / np.tan(
        np.radians(alt_deg + 7.31 / (alt_deg + 4.4))
    )
    refraction_arcmin *= (PRESS_MB / 1010) * (283 / (273 + TEMP_C))
    return refraction_arcmin / 60.0


def parse_solve_result(line):
    """Parse one nominal/actual coordinate result line."""
    parts = line.split()
    if len(parts) < 8 or parts[3] != "actual:":
        raise ValueError("expected '<name>: RA Dec actual: RA Dec DATE TIME'")
    try:
        nominal_ra = float(parts[1])
        nominal_dec = float(parts[2])
        actual_ra = float(parts[4])
        actual_dec = float(parts[5])
        observation_time = Time(f"{parts[6]} {parts[7]}")
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid solve result: {error}") from error
    return nominal_ra, nominal_dec, actual_ra, actual_dec, observation_time


def read_solve_results(filename):
    """Read all solve-result rows without discarding headerless observations."""
    path = Path(filename)
    records = []
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                records.append(parse_solve_result(line))
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: {error}") from error
    if not records:
        raise ValueError(f"{path}: no solve results found")
    return records


def wrapped_angle_difference(actual, nominal):
    """Return ``actual - nominal`` wrapped to [-pi, pi]."""
    difference = np.asarray(actual) - np.asarray(nominal)
    return np.arctan2(np.sin(difference), np.cos(difference))


def _coordinates_to_altaz(nominal_ra, nominal_dec, actual_ra, actual_dec, time):
    location = EarthLocation(lat=LAT * u.deg, lon=LON * u.deg)
    frame = AltAz(obstime=time, location=location)
    nominal = SkyCoord(
        ra=nominal_ra * u.deg,
        dec=nominal_dec * u.deg,
        frame="icrs",
    ).transform_to(frame)
    actual = SkyCoord(
        ra=actual_ra * u.deg,
        dec=actual_dec * u.deg,
        frame="icrs",
    ).transform_to(frame)
    return nominal, actual


def get_alt_az_offsets(line):
    """Convert one solve-result row to nominal Alt/Az and pointing offsets."""
    nominal_ra, nominal_dec, actual_ra, actual_dec, time = parse_solve_result(line)
    nominal, actual = _coordinates_to_altaz(
        nominal_ra,
        nominal_dec,
        actual_ra,
        actual_dec,
        time,
    )
    delta_alt = actual.alt.deg - nominal.alt.deg
    delta_az = np.degrees(
        wrapped_angle_difference(actual.az.rad, nominal.az.rad)
    ) * np.cos(nominal.alt.rad)
    return nominal.alt.deg, nominal.az.deg, delta_alt, delta_az


def result2dazel(filename):
    """Append calculated Alt/Az offsets and write beside the result input."""
    input_path = Path(filename)
    output_path = input_path.with_name(f"{input_path.stem}_dazel.txt")
    with input_path.open("r", encoding="utf-8") as source:
        with output_path.open("w", encoding="utf-8") as destination:
            for line_number, line in enumerate(source, start=1):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                try:
                    alt, az, delta_alt, delta_az = get_alt_az_offsets(line)
                except ValueError as error:
                    raise ValueError(
                        f"{input_path}:{line_number}: {error}"
                    ) from error
                destination.write(
                    line.strip()
                    + f"   {alt:7.4f} {az:8.4f} {delta_az:7.4f} {delta_alt:7.4f}\n"
                )
    return str(output_path)


def get_alt_az_data(file_path):
    """Load solve results as nominal Alt/Az and observed pointing offsets."""
    alts, azs, delta_alts, delta_azs = [], [], [], []
    for nominal_ra, nominal_dec, actual_ra, actual_dec, time in read_solve_results(
        file_path
    ):
        nominal, actual = _coordinates_to_altaz(
            nominal_ra,
            nominal_dec,
            actual_ra,
            actual_dec,
            time,
        )
        alts.append(np.radians(nominal.alt.deg - get_refraction_deg(nominal.alt.deg)))
        azs.append(nominal.az.rad)
        delta_alts.append(actual.alt.rad - nominal.alt.rad)
        delta_azs.append(
            wrapped_angle_difference(actual.az.rad, nominal.az.rad)
            * np.cos(nominal.alt.rad)
        )

    return tuple(
        np.asarray(values, dtype=float)
        for values in (alts, azs, delta_alts, delta_azs)
    )


def tpoint_residuals(coeffs, alts, azs, delta_alts, delta_azs):
    """Return weighted azimuth and elevation residuals for the TPOINT model."""
    ia, ie, ca, npae, an, aw, tf = coeffs
    model_delta_az = (
        -ia * np.cos(alts)
        - ca
        - npae * np.sin(alts)
        + an * np.sin(azs) * np.sin(alts)
        - aw * np.cos(azs) * np.sin(alts)
    )
    model_delta_alt = (
        ie + an * np.cos(azs) - aw * np.sin(azs) + tf * np.cos(alts)
    )
    residual_az = delta_azs - model_delta_az
    residual_alt = delta_alts - model_delta_alt
    return np.concatenate([residual_az, residual_alt])


def fit_pointing_model(alts, azs, delta_alts, delta_azs):
    """Fit seven TPOINT terms and validate that the solution is usable."""
    arrays = [
        np.asarray(values, dtype=float)
        for values in (alts, azs, delta_alts, delta_azs)
    ]
    lengths = {len(values) for values in arrays}
    if len(lengths) != 1:
        raise ValueError("pointing arrays must have equal lengths")
    if len(arrays[0]) < 4:
        raise ValueError("at least 4 observations are required for 7 coefficients")
    if not all(np.all(np.isfinite(values)) for values in arrays):
        raise ValueError("pointing arrays must contain only finite values")

    initial_guess = np.zeros(7)
    coefficients, covariance, _, message, status = leastsq(
        tpoint_residuals,
        initial_guess,
        args=tuple(arrays),
        full_output=True,
    )
    if status not in (1, 2, 3, 4):
        raise RuntimeError(f"pointing fit did not converge: {message}")
    if covariance is None:
        raise RuntimeError("pointing fit is rank-deficient")
    return coefficients


def radians_to_arcminutes(value):
    """Convert radians to arcminutes."""
    return np.degrees(value) * 60.0


def mountcal(filepath):
    """Fit, print, and plot the seven-term pointing calibration."""
    alts, azs, delta_alts, delta_azs = get_alt_az_data(filepath)
    coefficients = fit_pointing_model(alts, azs, delta_alts, delta_azs)
    labels = [
        "IA (Az Index)",
        "IE (El Index)",
        "CA (Collimation)",
        "NPAE (Axial Non-perp)",
        "AN (Az Tilt N/S)",
        "AW (Az Tilt E/W)",
        "TF (Tube Flexure)",
    ]

    print(f"{'Coefficient':<20} | {'Value (arcmin)':>15}")
    print("-" * 40)
    for label, value in zip(labels, coefficients):
        print(f"{label:<20} | {radians_to_arcminutes(value):>15.4f}")

    residuals = tpoint_residuals(
        coefficients,
        alts,
        azs,
        delta_alts,
        delta_azs,
    )
    observation_count = len(delta_azs)
    model_delta_azs = delta_azs - residuals[:observation_count]
    model_delta_alts = delta_alts - residuals[observation_count:]

    _, axes = plt.subplots(1, 2)
    axes[1].plot(np.degrees(alts), np.degrees(delta_alts), "o", label="observed")
    axes[1].plot(np.degrees(alts), np.degrees(model_delta_alts), "+", label="model")
    axes[1].set_xlabel("Altitude")
    axes[1].set_ylabel("Offset [deg]")
    axes[0].plot(np.degrees(azs), np.degrees(delta_azs), "o", label="observed")
    axes[0].plot(np.degrees(azs), np.degrees(model_delta_azs), "+", label="model")
    axes[0].set_xlabel("Azimuth")
    axes[0].set_ylabel("Offset [deg]")
    for axis in axes:
        axis.legend()
    return coefficients
