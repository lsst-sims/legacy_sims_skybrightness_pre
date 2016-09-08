import numpy as np

# Goal is to compute the approximate 5-sigma limiting depth at evenly spaced time intervals for each
# healpixel and field center. The we can have 6 histograms with ~200 bins, and then easily compute the
# approximate percentile of a spot in the sky at any time, making it easier to decide between
# fields or filters.

# Need to generate a histogram for every point in the sky...
__all__ = ['conditions2m5']


def conditions2m5(FWHMeff_zenith, sky_brightness, airmass, t_vis=30., filtername='r'):
    """
    Take the observing conditions and return the 5-sigma limiting depth

    Parameters
    ----------
    FWHMeff_zenith : float
        The seeing at zenith (arcsec)
    sky_brightness : float or np.array
        Sky background brightness (mag / sq arcsec)
    airmass : float or np.array
        The airmass(es) (unitless). Should be same length as sky_brightness
    t_vis : float (30.)
        The total exposure time of a visit (seconds)
    filtername : str
        The filter of the observation (ugrizy)

    Returns
    -------
    m5 : float or np.array
        The 5-sigma limiting depth of a point source taken in the conditions
    """

    # Only construct the dictionaries once
    if not hasattr(conditions2m5, 'C_m'):
        # Values from overview paper Table 2.
        # XXX--should check if there's a place where the C_m values can be computed from
        # all the latest throughput values
        conditions2m5.C_m = {'u': 22.92, 'g': 24.29, 'r': 24.33, 'i': 24.20, 'z': 24.07, 'y': 23.69}
        conditions2m5.k_m = {'u': 0.451, 'g': 0.163, 'r': 0.087, 'i': 0.065, 'z': 0.043, 'y': 0.138}

    # The FWHM at the airmass(es)
    FWHMeff = airmass**(0.6) * FWHMeff_zenith

    # Equation 6 from overview paper
    m5 = conditions2m5.C_m[filtername] + 0.5 * (sky_brightness - 21.) + 2.5 * np.log10(0.7 / FWHMeff)
    m5 += 1.25 * np.log10(t_vis / 30.) - conditions2m5.k_m[filtername] * (airmass - 1.)

    return m5
