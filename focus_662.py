import zwoasi as asi
import numpy as np
import matplotlib.pylab as plt

def init_camera(gain=252, exp_s=3):
    # Setup Camera
    env_filename = "/usr/local/lib/libASICamera2.so" 
    asi.init(env_filename)
    num_cameras = asi.get_num_cameras()
    if num_cameras == 0:
        raise ValueError("No cameras found")

    camera = asi.Camera(0)

    # Optimized Plate Solving Settings
    camera.set_control_value(asi.ASI_GAIN, gain)        # High sensitivity HCG mode
    exp_us = int(exp_s*1000000)
    camera.set_control_value(asi.ASI_EXPOSURE, exp_us)  # 3s is usually plenty for a 25mm lens
    camera.set_image_type(asi.ASI_IMG_RAW16)

    return camera

plt.ion()

def focus_image(exp_s=3, gain=252, scale='lin'):
    camera = init_camera(gain, exp_s)
    imnum = 1
    try:
        while True:
            buffer = camera.capture()
            image = np.frombuffer(buffer, dtype=np.uint16).reshape((1080, 1920))

            if imnum == 1:
                f, ax = plt.subplots()
                if scale == 'log':
                    im = ax.imshow(np.log10(image), cmap='gray')
                else:
                    im = ax.imshow(image, cmap='gray')
            else:
                if scale == 'log':
                    im.set_data(np.log10(image))
                else:
                    im.set_data(image)

            ax.set_title("Captured Image "+str(imnum))
            plt.draw()
            plt.pause(0.1)
            imnum += 1

    except KeyboardInterrupt:
        print('Close camera and exit')

    return image

def parse_args(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description='A script that snaps an exposure from a ZWO ASI120MM camera')
    parser.add_argument('--exp', type=float, help="Optional exposure time [s]")
    parser.add_argument('--gain', type=int, help="Optional gain [0-1000]")
    parser.add_argument('--scale', type=str, help="Use 'log' or 'lin' for scaling")
    return parser.parse_args(argv)

if __name__ == "__main__":
    args = parse_args()
    exp = 3
    if args.exp:
        exp = args.exp
    gain = 252
    if args.gain:
        gain = args.gain
    scale = 'lin'
    if args.scale:
        scale = args.scale

    print(focus_image(exp, gain=gain, scale=scale))
