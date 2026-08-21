import contextlib
import inspect
import io
import unittest
from unittest.mock import patch

from focus_662 import focus_image, parse_args


class Focus662CliBindingTest(unittest.TestCase):
    def test_cli_parses_gain_as_an_integer_for_the_camera_sdk(self):
        args = parse_args([
            "--exp", "0.01",
            "--gain", "252",
            "--scale", "log",
        ])

        self.assertEqual(args.exp, 0.01)
        self.assertEqual(args.gain, 252)
        self.assertIsInstance(args.gain, int)
        self.assertEqual(args.scale, "log")

    def test_cli_rejects_unknown_scale(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parse_args(["--scale", "sqrt"])

    def test_cli_arguments_bind_exposure_gain_and_scale(self):
        bound = inspect.signature(focus_image).bind(
            0.01,
            gain=252,
            scale="log",
        )

        self.assertEqual(
            bound.arguments,
            {"exp_s": 0.01, "gain": 252, "scale": "log"},
        )

    def test_camera_closes_when_capture_fails(self):
        class FailingCamera:
            def __init__(self):
                self.closed = False

            def capture(self):
                raise RuntimeError("capture failed")

            def close(self):
                self.closed = True

        camera = FailingCamera()
        with patch("focus_662.init_camera", return_value=camera):
            with self.assertRaisesRegex(RuntimeError, "capture failed"):
                focus_image(0.01, gain=252, scale="log")

        self.assertTrue(camera.closed)


if __name__ == "__main__":
    unittest.main()
