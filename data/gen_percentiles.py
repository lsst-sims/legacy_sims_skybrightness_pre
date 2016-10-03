import numpy as np
from lsst.sims.skybrightness_pre import SkyModelPre, conditions2m5
from lsst.sims.photUtils import LSSTdefaults


def generate_percentiles(nbins=20):
    """
    Make histograms of the 5-sigma limiting depths for each point and each filter.
    """

    sm = SkyModelPre()
    # Number of spatial pixels
    npix = sm.sb['r'].shape[-1]
    mid_mjds = (sm.mjd_left + sm.mjd_right) / 2.
    # One year should be fine.
    mid_mjds = [mid_mjds[0]]
    filters = ['u', 'g', 'r', 'i', 'z', 'y']

    histograms = np.zeros((nbins, npix), dtype=zip(filters, [float]*6))
    histogram_npts = np.zeros(npix, dtype=zip(filters, [int]*6))
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
                # Compute the mjds at the correct times
                m5s = np.interp(mjds, sm.info['mjds'], sm.sb[filtername][:, indx])
                m5s = m5s[np.isfinite(m5s)]
                m5s = np.sort(m5s)
                percentile_points = np.round(np.linspace(0, m5s.size-1, nbins))
                if m5s.size > percentile_points.size:
                    histograms[filtername][:, indx] = m5s[percentile_points.astype(int)]
                    histogram_npts[filtername][indx] = m5s.size
                # make the histogram for this point in the sky
                # histograms[filtername][:, indx] += np.histogram(m5s[np.isfinite(m5s)],
                #                                                bins=bins[filtername])[0]

    np.savez('percentile_m5_maps.npz', histograms=histograms, histogram_npts=histogram_npts)


if __name__ == '__main__':

    generate_percentiles()
