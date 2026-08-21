import inspect
import unittest

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


if __name__ == "__main__":
    unittest.main()
