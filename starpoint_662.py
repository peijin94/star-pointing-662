import zwoasi as asi
import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from time import sleep
import matplotlib.pylab as plt
from pathlib import Path

import shutil
import sys

def print_dual(left_text, right_text):
    # Get terminal width
    columns, _ = shutil.get_terminal_size()

    # Calculate how many spaces are needed between the two strings
    padding = columns - len(left_text) - len(right_text)

    # Ensure there's at least one space if the text is too long
    padding = max(padding, 1)

    # \r moves to start, then we print left + spaces + right
    sys.stdout.write(f"\r{left_text}{' ' * padding}{right_text}")
    sys.stdout.flush()

def startable2dict(filename):
    ''' Read star table file and convert list of stars to a list of "targets"
    '''
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    targets = []
    for line in lines[2:]:
        if len(line) == 0:
            break
        name = line[:4]
        utcstr = line.split()[-2]+' '+line.split()[-1]
        stime = Time(utcstr)
        rastr, decstr = line[17:41].split()
        ra_parts = [float(x) for x in rastr.split(':')]
        dec_parts = [float(x) for x in decstr.split(':')]
        is_negative = decstr[0].strip().startswith('-')
        targets.append({"name":name, "time": stime,
                        "ra_h": int(ra_parts[0]), "ra_m": int(ra_parts[1]),
                        "ra_s": ra_parts[2], "dec_neg": is_negative, "dec_d": int(dec_parts[0]), 
                        "dec_m": int(dec_parts[1]), "dec_s": dec_parts[2]
                       })
    return targets

def init_camera():
    # Setup Camera
    env_filename = "/usr/local/lib/libASICamera2.so"
    asi.init(env_filename)
    num_cameras = asi.get_num_cameras()
    if num_cameras == 0:
        raise ValueError("No cameras found")
    return asi.Camera(0)

def snap(srcname='test', ra=None, dec=None, gain=252, exp_s=3, show=True):
    # Capture and save one image from ZWO ASI662MC camera
    camera = init_camera()
    # Optimized Plate Solving Settings
    camera.set_control_value(asi.ASI_GAIN, gain)        # High sensitivity HCG mode
    exp_us = int(exp_s*1000000)
    camera.set_control_value(asi.ASI_EXPOSURE, exp_us)  # 3s is usually plenty for a 25mm lens
    camera.set_image_type(asi.ASI_IMG_RAW16)

    # Capture image buffer
    print('Exposing',exp_s,'seconds on',srcname)
    buffer = camera.capture()
    image_data = np.frombuffer(buffer, dtype=np.uint16).reshape((1080, 1920))

    # Create FITS with ASTAP-Specific Headers
    hdu = fits.PrimaryHDU(image_data)
    hdr = hdu.header

    # --- CRITICAL HEADERS FOR ASTAP ---
    hdr['FOCALLEN'] = 25.0          # Your lens focal length in mm
    if ra:
        hdr['RA']  = ra    # Approx RA (if known) speeds up solving
    if dec:
        hdr['DEC'] = dec   # Approx Dec (if known)
    hdr['INSTRUME'] = 'ASI662MC'
    hdr['GAIN']     = 252
    hdr['EXPTIME']  = exp_s
    hdr['XPIXSZ']   = 2.9           # Pixel size in microns (essential for FOV calc)
    hdr['YPIXSZ']   = 2.9

    if ra and dec:
        filename = srcname+'_'+str(int(ra))+'_'+str(int(dec))+'_'+str(exp_s)+'.fits'
    else:
        filename = srcname+'_unk_unk_'+str(int(exp_s))+'.fits'
    # Save file
    print('Writing file',filename)
    hdu.writeto(filename, overwrite=True)

    camera.close()
    if show:
        plt.ion()
        plt.imshow(np.log10(image_data),origin='lower',cmap='gray')
        plt.pause(1)

    return filename

import subprocess
import os

def plate_solve(filename):
    # Command to solve: 
    # -f: input file
    # -ra / -dec: approximate search center (optional but speeds it up)
    # -fov: field of view (7.2 deg for 25 mm lens)
    # -r: search radius (3 degrees for 25 mm lens)
    # -v: silent mode/version
    command = ["astap", "-f", filename, "-fov", "7.2", "-r", "3"]

    #print(f"Running ASTAP on {filename}...")
    result = subprocess.run(command, capture_output=True, text=True)

    # ASTAP creates a .wcs file with the same name as the input if it succeeds
    wcs_file = filename.replace(".fits", ".wcs")

    if os.path.exists(wcs_file):
        #print("Solve Successful!")
        with open(wcs_file, 'r') as f:
            lines = f.readlines()
            # Extract basic coordinates from the WCS header file
            solution = {line.split('=')[0].strip(): line.split('=')[1].split('/')[0].strip() 
                        for line in lines if '=' in line}

            ra = float(solution.get('CRVAL1'))
            dec = float(solution.get('CRVAL2'))
            #print(f"Center RA: {ra}, Dec: {dec}")
            return ra, dec
    else:
        #print("Solve Failed. Check your FOV and star database.")
        return None

def imaging(startable, exp_s=3, test=False):
    ''' This does everything, from reading the startable file, taking an image,
        and doing the plate solve.
    '''
    targets = startable2dict(startable)
    dt = TimeDelta(30 * u.s)  # 30 s (into the future)

    datstr, timstr = targets[0]['time'].iso.split()
    timstr = timstr.replace(':','')

    # Make a new folder with the date of the first target
    folder = Path(datstr)
    if not folder.is_dir():
        os.mkdir(folder)
    os.chdir(folder)

    with open(datstr+'_solve_results.txt','w') as f:
        for target in targets:
            # This line is in the future, so wait next time and report out
            while Time.now() < target['time'] + dt:
                sleep(1)
                print_dual('Waiting for '+target['name']+' time '+target['time'].iso[:19],Time.now().iso[:19])
            if Time.now() >= target['time'] + dt:
                # We should have been moving to the target for at least 30 s
                # Check whether to do this line or skip it
                if Time.now() < target['time'] + 2*dt:
                    # We are in the 30-s window, so snap an image and plate solve it
                    ra = target['ra_h']*15 + target['ra_m']/4. + target['ra_s']/240.
                    if target['dec_neg']:
                        dec = target['dec_d'] - target['dec_m']/60. - target['dec_s']/3600.
                    else:
                        dec = target['dec_d'] + target['dec_m']/60. + target['dec_s']/3600.
                    if test:
                        print(f"{target['name']}: {ra:8.4f} {dec:8.4f} test")
                    else:
                        timstr = target['time'].iso.split()[1]
                        timstr = timstr.replace(':','')
                        fitsfile = snap(srcname=target['name']+'_'+timstr[:4], ra=ra, dec=dec, gain=252, exp_s=exp_s, show=False)
                        res = plate_solve(fitsfile)
                        if res:
                            ra_c, dec_c = res
                            success_line = f"{target['name']}: {ra:8.4f} {dec:8.4f} actual: {ra_c:8.4f} {dec_c:8.4f} {target['time'].iso[:19]}"
                            print(success_line)
                            f.write(success_line+'\n')
                            f.flush()
                        else:
                            print(target['name'],'solve failed!')
                    # Sleep until next line:
                    while Time.now() < target['time'] + dt*2:
                        sleep(1)
                        print_dual('Completed '+target['name']+' time '+target['time'].iso[:19],Time.now().iso[:19])
                    print()
                else:
                    # This line is more than 1 min into the past, so skip it and go to the next target
                    print('Skipped',target['name'],target['time'].iso[:19])

    print('All Done!')

from astropy.time import Time
import astropy.units as u

# --- CONFIGURATION ---
LAT = 37.5177  # Your latitude
LON = -118.3655 # Your longitude (negative for West)
TEMP_C = 10.0   # Ground temperature (Celsius)
PRESS_MB = 870 # Pressure (millibars)
# ---------------------

def get_refraction_deg(alt_deg):
    """
    Bennett's formula for atmospheric refraction in degrees.
    """
    if alt_deg <= 0: return 0
    # Refraction in arcminutes
    R = 1.0/np.tan(np.radians(alt_deg + (7.31/(alt_deg + 4.4))))
    
    # Correct for Temp/Pressure (approximate)
    R *= (PRESS_MB/1010) * (283/(273 + TEMP_C))
    return R/60.0

def get_alt_az_offsets(line):
    # Example line parsing
    parts = line.split()
    ra_nom, dec_nom = float(parts[1]), float(parts[2])
    ra_act, dec_act = float(parts[4]), float(parts[5])
    obs_time = f"{parts[6]} {parts[7]}"

    location = EarthLocation(lat=LAT*u.deg, lon=LON*u.deg)
    time = Time(obs_time)

    # Define the two points in RA/Dec
    coord_nom = SkyCoord(ra=ra_nom*u.deg, dec=dec_nom*u.deg, frame='icrs')
    coord_act = SkyCoord(ra=ra_act*u.deg, dec=dec_act*u.deg, frame='icrs')

    # Transform both to AltAz for that specific time and location
    altaz_frame = AltAz(obstime=time, location=location)
    azalt_nom = coord_nom.transform_to(altaz_frame)
    azalt_act = coord_act.transform_to(altaz_frame)

    # Calculate differences
    d_alt = azalt_act.alt.deg - azalt_nom.alt.deg
    # Correct Azimuth for the cosine of Elevation (similar to RA/Dec)
    d_az = (azalt_act.az.deg - azalt_nom.az.deg) * np.cos(np.radians(azalt_nom.alt.deg))

    return azalt_nom.alt.deg, azalt_nom.az.deg, d_alt, d_az

def result2dazel(filename):
    basename = os.path.basename(filename)
    with open(basename[:-4]+'_dazel.txt','w') as o:
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                alt, az, delta_alt, delta_az = get_alt_az_offsets(line)
                line = line.strip()+f'   {alt:7.4f} {az:8.4f} {delta_az:7.4f} {delta_alt:7.4f}'
                o.write(line+'\n')

from scipy.optimize import leastsq

# --- CONFIGURATION ---
FILE_PATH = "2026-04-24_results.txt"

def get_alt_az_data(file_path):
    alts, azs, d_alts, d_azs = [], [], [], []
    location = EarthLocation(lat=LAT*u.deg, lon=LON*u.deg)

    with open(file_path, 'r') as f:
        next(f) # Skip header
        for line in f:
            if not line.strip(): continue
            parts = line.split()
            # Parsing: idx, ra_nom, dec_nom, 'actual:', ra_act, dec_act, date, time
            ra_nom, dec_nom = float(parts[1]), float(parts[2])
            ra_act, dec_act = float(parts[4]), float(parts[5])
            t = Time(f"{parts[6]} {parts[7]}")

            # Convert to Alt/Az
            frame = AltAz(obstime=t, location=location)
            nom = SkyCoord(ra=ra_nom*u.deg, dec=dec_nom*u.deg, frame='icrs').transform_to(frame)
            act = SkyCoord(ra=ra_act*u.deg, dec=dec_act*u.deg, frame='icrs').transform_to(frame)

            alts.append(np.radians(nom.alt.deg - get_refraction_deg(nom.alt.deg)))
            azs.append(nom.az.rad)
            # Observed minus Calculated
            d_alts.append(act.alt.rad - nom.alt.rad)
            # Azimuth difference must be weighted by cos(elevation)
            d_azs.append((act.az.rad - nom.az.rad) * np.cos(nom.alt.rad))

    return np.array(alts), np.array(azs), np.array(d_alts), np.array(d_azs)

def tpoint_residuals(coeffs, alts, azs, d_alts, d_azs):
    IA, IE, CA, NPAE, AN, AW, TF = coeffs
    
    # TPOINT geometric model equations for Alt-Az
    # Predicted Azimuth error (weighted by cos(alt))
    mod_d_az = (-IA * np.cos(alts) + 
                -CA + 
                -NPAE * np.sin(alts) + 
                AN * np.sin(azs) * np.sin(alts) - 
                AW * np.cos(azs) * np.sin(alts))

    # Predicted Elevation error
    mod_d_alt = (IE + 
                 AN * np.cos(azs) + 
                 -AW * np.sin(azs) + 
                 TF * np.cos(alts))

    res_az = d_azs - mod_d_az
    res_alt = d_alts - mod_d_alt
    return np.concatenate([res_az, res_alt])

def mountcal(filepath):
    # Load and Fit
    alts, azs, d_alts, d_azs = get_alt_az_data(filepath)
    initial_guess = [0, 0, 0, 0, 0, 0, 0] # Start at zero for all terms
    final_coeffs, success = leastsq(tpoint_residuals, initial_guess, args=(alts, azs, d_alts, d_azs))

    # --- OUTPUT RESULTS ---
    labels = ["IA (Az Index)", "IE (El Index)", "CA (Collimation)", 
              "NPAE (Axial Non-perp)", "AN (Az Tilt N/S)", "AW (Az Tilt E/W)", "TF (Tube Flexure)"]

    print(f"{'Coefficient':<20} | {'Value (arcmin)':>15}")
    print("-" * 40)
    for label, val in zip(labels, final_coeffs):
        # Convert radians to degrees for readability
        val_arcmin = np.degrees(val)
        print(f"{label:<20} | {val_arcmin:>15.4f}")
       
    result = tpoint_residuals(final_coeffs, alts, azs, d_alts, d_azs)
    f, ax = plt.subplots(1,2)
    nazs = len(d_azs)
    ax[1].plot(np.degrees(alts), np.degrees(d_alts), 'o')
    ax[1].plot(np.degrees(alts), np.degrees(result[nazs:]), '+')
    ax[1].set_xlabel('Altitude')
    ax[1].set_ylabel('Offset [deg]')
    ax[0].plot(np.degrees(azs), np.degrees(d_azs), 'o')
    ax[0].plot(np.degrees(azs), np.degrees(result[:nazs]), '+')
    ax[0].set_xlabel('Azimuth')
    ax[0].set_ylabel('Offset [deg]')
    return final_coeffs
