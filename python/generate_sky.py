import numpy as np
import lsst.sims.skybrightness as sb
import lsst.sims.utils as utils
import healpy as hp
import os, sys

# XXX--switch this up to dynamically decide what the timestep should be.  Say, if linear interpolation gets the right answer
# to within 0.3 mags, then it's good enough? That way we get smaller timesteps at moon rise/set and twilight, longer when stable!


def generate_sky(mjd0=59560.2, duration=0.01, timestep=5., outfile='generated_sky.npz', nside=32,
                 sunLimit=-12., fieldID=False, airmass_limit=1.5, dm=0.3, verbose=True):
    """
    Use the sky brightness model to generate a number of useful numpy arrays that can be used
    to look-up sky brighntess and other pre-computed info
    """

    sunLimit = np.radians(sunLimit)

    # Set the time steps
    mjds = []
    mjd_max = mjd0 + duration
    timestep = timestep / 60. / 24.  # Convert to days

    if verbose:
        print 'generating %i unique MJDs' % mjds.size

    # Maybe add in relevant time points like moon rise/set and twilight times

    # set the RA, Dec array
    hpindx = np.arange(hp.nside2npix(nside))
    ra, dec = utils.hpid2RaDec(nside, hpindx)

    # Switch the indexing to opsim field ID if requested
    if fieldID:
        field_data = np.loadtxt('fieldID.dat', delimiter='|', skiprows=1,
                                dtype=zip(['id', 'ra', 'dec'],[int, float, float]))
        ra = np.radians(field_data['ra'])
        dec = np.radians(field_data['dec'])
    else:
        hpindx = np.arange(hp.nside2npix(nside))
        ra, dec = utils.hpid2RaDec(nside, hpindx)

    if verbose:
        print 'using %i points on the sky' % ra.size

    # Set up the sky brightness model
    sm = sb.SkyModel(mags=True)

    sky_properties = {}
    filter_names = ['u', 'g', 'r', 'i', 'z', 'y']

    names = ['sky_brightness'] #, 'nominal_m5']
    types = [float] * len(names)

    for filter_name in filter_names:
        sky_properties[filter_name] = np.zeros((mjds.size, ra.size), dtype=float)  #  dtype=zip(names, types))

    names = ['airmass'] #, 'd2moon', 'mask']
    types = [float] #, float, int]
    general_properties = np.zeros((mjds.size, ra.size))  #  , zip(names, types))
    sunAlts = np.zeros(mjds.size, dtype=float)
    sunAlts.fill(-666)

    maxI = mjds.size
    for i, mjd in enumerate(mjds):
        sm.setRaDecMjd(ra, dec, mjd, degrees=True)
        sunAlts[i] = sm.sunAlt
        if sm.sunAlt <= sunLimit:
            mags = sm.returnMags()
            for key in filter_names:
                sky_properties[key][i, :] += mags[key]
            general_properties[i, :] += sm.airmass
        progress = i/float(maxI)*100
        text = "\rprogress = %.1f%%"%progress
        sys.stdout.write(text)
        sys.stdout.flush()
            # XXX-need to call a haversine here I think
            #general_properties['d2moon'][:, i] =
        print ''


    toSave = np.where(sunAlts <= np.radians(sunLimit))[0]
    # Collapse down and remove the mjds that had no info
    for key in filter_names:
        sky_properties[key][np.isnan(sky_properties[key])] = 0
        sky_properties[key] = sky_properties[key][toSave, :]
        sky_properties[key][np.where(sky_properties[key] == 0)] = hp.UNSEEN

    mjds = mjds[toSave]
    sunAlts = sunAlts[toSave]
    # XXX--collapse the arrays down and only save those with data
    maxes = np.zeros((toSave.size, len(filter_names)), dtype=float)
    for i, filter_name in enumerate(filter_names):
        maxes[:, i] = np.max(sky_properties[:, i], dimension=1)
    maxes = np.max(maxes, dimension=1)

    toSave = np.where(maxes > 0)[0]
    for key in filter_names:
        sky_properties[key] = sky_properties[key][toSave, :]
    mjds = mjds[toSave]
    sunAlts = sunAlts[toSave]



    np.savez(outfile, general_properties=general_properties,
             sky_properties=sky_properties, mjds=mjds, sunAlts=sunAlts)

if __name__ == "__main__":
    generate_sky()

