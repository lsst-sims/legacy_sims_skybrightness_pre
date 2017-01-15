import numpy as np
from calcM5 import calcM5
import bandpassUtils as bu
from lsst.utils import getPackageDir

# Need to setup syseng_throughputs

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
        defaultDirs = bu.setDefaultDirs(rootDir = getPackageDir('syseng_throughputs'))
        addLosses = True
        atmosphere = bu.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
        hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses, atmosphereOverride=atmosphere)
        m5 = calcM5(hardware, system, atmosphere, title='', return_t2_values=True)
        conditions2m5.C_m = m5['Cm']  # {'u': 22.92, 'g': 24.29, 'r': 24.33, 'i': 24.20, 'z': 24.07, 'y': 23.69}
        conditions2m5.dCm_inf = m5['dCm_infinity']  # {'u': 0.67, 'g': 0.21, 'r': 0.11, 'i': 0.08, 'z': 0.05, 'y': 0.04}
        conditions2m5.k_m = m5['kAtm']  # {'u': 0.451, 'g': 0.163, 'r': 0.087, 'i': 0.065, 'z': 0.043, 'y': 0.138}

    # The FWHM at the airmass(es)
    FWHMeff = airmass**(0.6) * FWHMeff_zenith

    # Equation 6 from overview paper
    m5 = conditions2m5.C_m[filtername] + 0.5 * (sky_brightness - 21.) + 2.5 * np.log10(0.7 / FWHMeff)
    m5 += 1.25 * np.log10(t_vis / 30.) - conditions2m5.k_m[filtername] * (airmass - 1.)

    # Equation 7 from the overview paper
    if t_vis > 30.:
        tau = t_vis/30.
        numerator = 10.**(0.8 * conditions2m5.dCm_inf[filtername]) - 1.
        dcm = conditions2m5.dCm_inf[filtername]-1.25*np.log10(1. + numerator/tau)
        m5 += dcm

    return m5
