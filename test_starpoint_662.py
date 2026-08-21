import contextlib
import io
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.utils import iers

import starpoint_662 as starpoint


iers.conf.auto_download = False


class FakeCamera:
    def __init__(self, buffer=None, error=None):
        self.buffer = buffer
        self.error = error
        self.closed = False
        self.controls = []
        self.image_type = None

    def set_control_value(self, control, value):
        self.controls.append((control, value))

    def set_image_type(self, image_type):
        self.image_type = image_type

    def capture(self):
        if self.error is not None:
            raise self.error
        return self.buffer

    def close(self):
        self.closed = True


class StarTableTest(unittest.TestCase):
    def test_valid_table_parses_negative_subdegree_declination(self):
        row = (
            "STAR"
            + " " * 13
            + "12:34:56.0 -00:30:00.0  scheduled 2026-08-20 12:00:00\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory, "targets.txt")
            table.write_text("header one\nheader two\n" + row, encoding="utf-8")

            target = starpoint.startable2dict(table)[0]

        self.assertEqual(target["name"], "STAR")
        self.assertAlmostEqual(
            starpoint.target_coordinates_deg(target)[0],
            188.73333333333332,
        )
        self.assertAlmostEqual(starpoint.target_coordinates_deg(target)[1], -0.5)

    def test_empty_table_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory, "empty.txt")
            table.write_text("header one\nheader two\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "no targets"):
                starpoint.startable2dict(table)

    def test_malformed_row_reports_its_line(self):
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory, "bad.txt")
            table.write_text("header one\nheader two\nnot a target\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, r"bad\.txt:3"):
                starpoint.startable2dict(table)


class CameraCaptureTest(unittest.TestCase):
    def test_zero_coordinates_and_requested_gain_reach_fits_header(self):
        image = np.zeros((1080, 1920), dtype=np.uint16)
        camera = FakeCamera(buffer=image.tobytes())

        with tempfile.TemporaryDirectory() as directory:
            source = str(Path(directory, "zero"))
            with patch("starpoint_662.init_camera", return_value=camera):
                filename = starpoint.snap(
                    source,
                    ra=0.0,
                    dec=0.0,
                    gain=123,
                    exp_s=0.01,
                    show=False,
                )

            header = fits.getheader(filename)
            self.assertEqual(header["RA"], 0.0)
            self.assertEqual(header["DEC"], 0.0)
            self.assertEqual(header["GAIN"], 123)
            self.assertNotIn("_unk_unk_", filename)
            self.assertTrue(camera.closed)

    def test_camera_closes_when_capture_fails(self):
        camera = FakeCamera(error=RuntimeError("capture failed"))

        with patch("starpoint_662.init_camera", return_value=camera):
            with self.assertRaisesRegex(RuntimeError, "capture failed"):
                starpoint.snap(show=False)

        self.assertTrue(camera.closed)


class PlateSolveTest(unittest.TestCase):
    def test_failed_solver_cannot_reuse_stale_wcs(self):
        with tempfile.TemporaryDirectory() as directory:
            fits_path = Path(directory, "image.fits")
            wcs_path = fits_path.with_suffix(".wcs")
            fits_path.touch()
            wcs_path.write_text("CRVAL1=12.5\nCRVAL2=-3.5\n", encoding="utf-8")

            result = subprocess.CompletedProcess([], 1, "", "failed")
            with patch("starpoint_662.subprocess.run", return_value=result):
                self.assertIsNone(starpoint.plate_solve(fits_path, timeout_s=5))

            self.assertFalse(wcs_path.exists())

    def test_successful_solver_coordinates_are_parsed(self):
        with tempfile.TemporaryDirectory() as directory:
            fits_path = Path(directory, "image.fits")
            fits_path.touch()

            def solve(*args, **kwargs):
                fits_path.with_suffix(".wcs").write_text(
                    "CRVAL1 = 123.25 / degrees\nCRVAL2 = -45.5 / degrees\n",
                    encoding="utf-8",
                )
                return subprocess.CompletedProcess(args, 0, "", "")

            with patch("starpoint_662.subprocess.run", side_effect=solve) as run:
                self.assertEqual(
                    starpoint.plate_solve(fits_path, timeout_s=7),
                    (123.25, -45.5),
                )

            self.assertEqual(run.call_args.kwargs["timeout"], 7)


class ImagingPathTest(unittest.TestCase):
    def test_imaging_does_not_change_the_process_directory(self):
        target = {
            "name": "TEST",
            "time": Time.now() - TimeDelta(120 * u.s),
            "ra_h": 1,
            "ra_m": 2,
            "ra_s": 3.0,
            "dec_neg": False,
            "dec_d": 4,
            "dec_m": 5,
            "dec_s": 6.0,
        }
        original = Path.cwd()

        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            try:
                before = Path.cwd()
                with patch("starpoint_662.startable2dict", return_value=[target]):
                    starpoint.imaging("unused.txt", test=True)
                self.assertEqual(Path.cwd(), before)
            finally:
                os.chdir(original)


class PointingCalibrationTest(unittest.TestCase):
    def test_headerless_results_keep_the_first_observation(self):
        rows = [
            "A001: 10.0 20.0 actual: 10.1 20.1 2026-08-20 12:00:00",
            "A002: 30.0 40.0 actual: 30.1 40.1 2026-08-20 12:05:00",
        ]
        with tempfile.TemporaryDirectory() as directory:
            result_file = Path(directory, "results.txt")
            result_file.write_text("\n".join(rows) + "\n", encoding="utf-8")

            records = starpoint.read_solve_results(result_file)

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0][0], 10.0)

    def test_azimuth_delta_wraps_across_zero(self):
        actual = np.radians(0.2)
        nominal = np.radians(359.8)

        wrapped = starpoint.wrapped_angle_difference(actual, nominal)

        self.assertAlmostEqual(np.degrees(wrapped), 0.4, places=10)

    def test_radians_are_converted_to_arcminutes(self):
        self.assertAlmostEqual(
            starpoint.radians_to_arcminutes(np.radians(1.0)),
            60.0,
        )

    def test_fit_requires_enough_observations(self):
        values = np.zeros(3)
        with self.assertRaisesRegex(ValueError, "at least 4"):
            starpoint.fit_pointing_model(values, values, values, values)

    def test_fit_recovers_synthetic_coefficients(self):
        alts = np.radians(np.linspace(15, 75, 20))
        azs = np.radians(np.linspace(3, 350, 20))
        expected = np.array([1, -2, 3, -4, 5, -6, 7], dtype=float) * 1e-4
        model = -starpoint.tpoint_residuals(
            expected,
            alts,
            azs,
            np.zeros_like(alts),
            np.zeros_like(alts),
        )
        split = len(alts)

        actual = starpoint.fit_pointing_model(
            alts,
            azs,
            model[split:],
            model[:split],
        )

        np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=1e-10)

    def test_dazel_output_is_written_beside_input(self):
        row = "A001: 10.0 20.0 actual: 10.1 20.1 2026-08-20 12:00:00\n"
        with tempfile.TemporaryDirectory() as directory:
            result_file = Path(directory, "results.txt")
            result_file.write_text(row, encoding="utf-8")

            with patch(
                "starpoint_662.get_alt_az_offsets",
                return_value=(30.0, 40.0, 0.1, 0.2),
            ):
                output = starpoint.result2dazel(result_file)

            self.assertEqual(Path(output), Path(directory, "results_dazel.txt"))
            self.assertTrue(Path(output).exists())


if __name__ == "__main__":
    with contextlib.redirect_stderr(io.StringIO()):
        unittest.main()
