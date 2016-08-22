import numpy as np
import lsst.sims.skybrightness as sb
import lsst.sims.utils as utils
import healpy as hp
import sys
import ephem
from lsst.sims.skybrightness.utils import mjd2djd


def generate_sky(mjd0=59560.2, duration=0.01, timestep=5., timestep_max=20.,
                 outfile='generated_sky.npz', nside=32,
                 sunLimit=-12., fieldID=False, airmass_limit=1.5, dm=0.3, verbose=True):
    """
    Use the sky brightness model to generate a number of useful numpy arrays that can be used
    to look-up sky brighntess and other pre-computed info
    """

    sunLimit = np.radians(sunLimit)

    # Set the time steps
    mjd_max = mjd0 + duration*365.25
    timestep = timestep / 60. / 24.  # Convert to days
    timestep_max = timestep_max / 60. / 24.  # Convert to days
    # Switch the indexing to opsim field ID if requested

    # Look at the mjds and toss ones where the sun is up
    mjds = np.arange(mjd0, mjd_max+timestep, timestep)
    sunAlts = np.zeros(mjds.size, dtype=float)

    telescope = utils.Site('LSST')
    Observatory = ephem.Observer()
    Observatory.lat = telescope.latitude_rad
    Observatory.lon = telescope.longitude_rad
    Observatory.elevation = telescope.height

    sun = ephem.Sun()

    # Compute the sun altitude for all the possible MJDs
    for i, mjd in enumerate(mjds):
        Observatory.date = mjd2djd(mjd)
        sun.compute(Observatory)
        sunAlts[i] = sun.alt

    mjds = mjds[np.where(sunAlts <= np.radians(sunLimit))]


    if fieldID:
        field_data = np.loadtxt('fieldID.dat', delimiter='|', skiprows=1,
                                dtype=zip(['id', 'ra', 'dec'], [int, float, float]))
        ra = field_data['ra']
        dec = field_data['dec']
    else:
        hpindx = np.arange(hp.nside2npix(nside))
        ra, dec = utils.hpid2RaDec(nside, hpindx)

    if verbose:
        print 'using %i points on the sky' % ra.size
        print 'using %i mjds' % mjds.size

    # Set up the sky brightness model
    sm = sb.SkyModel(mags=True)

    filter_names = ['u', 'g', 'r', 'i', 'z', 'y']

    # Initialize the relevant lists
    dict_of_lists = {'airmass': [], 'sunAlts': [], 'mjds': []}
    sky_brightness = {}
    for filter_name in filter_names:
        sky_brightness[filter_name] = []

    length = mjds[-1] - mjds[0]

    for mjd in mjds:
        progress = (mjd-mjd0)/length*100
        text = "\rprogress = %.1f%%"%progress
        sys.stdout.write(text)
        sys.stdout.flush()
        sm.setRaDecMjd(ra, dec, mjd, degrees=True)
        if sm.sunAlt <= sunLimit:
            mags = sm.returnMags()
            for key in filter_names:
                sky_brightness[key].append(mags[key])
            dict_of_lists['airmass'].append(sm.airmass)
            dict_of_lists['sunAlts'].append(sm.sunAlt)
            dict_of_lists['mjds'].append(mjd)
            # XXX don't toss things now.
            if False:
                # Check if we can interpolate the second to last sky brightnesses
                overhead = np.where((dict_of_lists['airmass'][-1] <= airmass_limit) &
                                    (dict_of_lists['airmass'][-2] <= airmass_limit))
                # XXX -- might need to also excude things near the moon
                if np.size(overhead[0]) > 0:
                    can_interp = True
                    mjd1 = dict_of_lists['mjds'][-3]
                    mjd2 = dict_of_lists['mjds'][-2]
                    mjd3 = dict_of_lists['mjds'][-1]
                    # Linear interpolation weights
                    wterm = (mjd2 - mjd1) / (mjd3-mjd1)
                    w1 = 1. - wterm
                    w2 = wterm
                    for filter_name in filter_names:
                        interp_sky = w1 * sky_brightness[filter_name][-3][overhead]
                        interp_sky += w2 * sky_brightness[filter_name][-1][overhead]
                        diff = np.abs(sky_brightness[filter_name][-2][overhead]-interp_sky)
                        if np.size(diff[~np.isnan(diff)]) > 0:
                            if np.max(diff[~np.isnan(diff)]) > dm:
                                can_interp = False
                    # If the timestep is getting big, don't interp even if you can
                    if (dict_of_lists['mjds'][-3] - dict_of_lists['mjds'][-2]) > timestep_max:
                        can_interp = False
                    if can_interp:
                        for key in dict_of_lists:
                            del dict_of_lists[key][-2]
                        for key in sky_brightness:
                            del sky_brightness[key][-2]
    print ''

    for key in dict_of_lists:
        dict_of_lists[key] = np.array(dict_of_lists[key])
    for key in sky_brightness:
        sky_brightness[key] = np.array(sky_brightness[key])

    np.savez(outfile, dict_of_lists = dict_of_lists, sky_brightness=sky_brightness)

if __name__ == "__main__":

    #generate_sky(fieldID=True, outfile='generated_sky_field.npz')
    generate_sky()
