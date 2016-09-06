import numpy as np
import LSSTdefaults

# Goal is to compute the approximate 5-sigma limiting depth at evenly spaced time intervals for each
# healpixel and field center. The we can have 6 histograms with ~200 bins, and then easily compute the
# approximate percentile of a spot in the sky at any time, making it easier to decide between
# fields or filters.

# Need to generate a histogram for every point in the sky...

def generate_percentile(FWHMeff=0.80, filtername='r'):
    """
    Load all the available dates and compute the percentile conditions for each healpixel and field

    Parameters
    ----------
    fwhm_effective : float
        The zenith fwhm effective value to assume for computing 5-sigma values
    """

    # Compute the 5-sigma limiting depth maps.

    # 

    if FWHMeff is None:
        FWHMeff = LSSTdefaults().FWHMeff(filtername)