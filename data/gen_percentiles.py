import numpy as np
from lsst.sims.skybrightness_pre import SkyModelPre, conditions2m5
from lsst.sims.photUtils import LSSTdefaults


def generate_percentiles(nbins=200):
    """
    Make histograms of the 5-sigma limiting depths for each point and each filter.
    """

    sm = SkyModelPre()
    # Number of spatial pixels
    npix = sm.sb['r'].shape[-1]
    mid_mjds = (sm.mjd_left + sm.mjd_right) / 2.
    filters = ['u', 'g', 'r', 'i', 'z', 'y']

    bins = np.zeros(nbins+1, dtype=zip(filters, [float]*6))
    mag_min = 16
    mag_max = 25
    for filtername in filters:
        bins[filtername] = np.linspace(mag_min, mag_max, num=nbins+1)

    histograms = np.zeros((nbins, npix), dtype=zip(filters, [int]*6))

    timestep = 10.  # minutes
    timestep = timestep / 60. / 24.

    # Should pull this value from the saved header info in the future
    ts_max = 20.
    ts_max = ts_max / 60. / 24.

    for mjd in mid_mjds:
        sm._load_data(mjd)
        # make an array of mjds evenly spaced
        all_mjds = np.arange(sm.info['mjds'].min(), sm.info['mjds'].max(), timestep)
        left = np.searchsorted(sm.info['mjds'], all_mjds)
        right = np.searchsorted(sm.info['mjds'], all_mjds, side='right')
        d1 = sm.info['mjds'][left] - all_mjds
        d2 = sm.info['mjds'][right] - all_mjds

        # mask out the mjds if the nearest point is masked.
        mjd_indxes = left.copy()
        other_closer = np.where(d2 < d1)
        mjd_indxes[other_closer] = right[other_closer]
        gaps = sm.info['mjds'][right] - sm.info['mjds'][right-1]
        mjds = all_mjds[np.where(gaps < ts_max)]

        for filtername in filters:
            # convert surface brightness to m5
            FWHMeff = LSSTdefaults().FWHMeff(filtername)
            sm.sb[filtername] = conditions2m5(FWHMeff, sm.sb[filtername],
                                              sm.info['airmass'], filtername=filtername)
            for indx in np.arange(npix):
                # Toss out mjds that should be masked
                # mjds = all_mjds[sm.info['masks'][mjd_indxes, indx]]
                m5s = np.interp(mjds, sm.info['mjds'], sm.sb[filtername][:, indx])
                # make the histogram for this point in the sky
                histograms[filtername][:, indx] += np.histogram(m5s[np.isfinite(m5s)],
                                                                bins=bins[filtername])[0]

    np.savez('percentile_m5_maps.npz', histograms=histograms, bins=bins)


if __name__ == '__main__':

    generate_percentiles()
