import unittest
import lsst.sims.skybrightness as sb
import lsst.sims.skybrightness_pre as sbp
import lsst.utils.tests
import healpy as hp
import numpy as np
import lsst.sims.utils as utils
import ephem
from lsst.sims.skybrightness.utils import mjd2djd
from lsst.sims.utils import haversine


class TestSkyPre(unittest.TestCase):

    def testSunMoon(self):
        """
        Test that the sun moon interpolation is good enough
        """
        sm = sbp.SkyModelPre()
        telescope = utils.Site('LSST')
        Observatory = ephem.Observer()
        Observatory.lat = telescope.latitude_rad
        Observatory.lon = telescope.longitude_rad
        Observatory.elevation = telescope.height

        sun = ephem.Sun()
        moon = ephem.Moon()

        mjd1 = sm.info['mjds'][0]
        mjd2 = sm.info['mjds'][3]

        mjds = np.linspace(mjd1, mjd2, 20)

        # Demand Moon and Sun Positions match to within 3 arcmin
        arcmin_places = np.abs(np.floor(np.log10(3./60./180.*np.pi))).astype(int)

        for mjd in mjds:
            Observatory.date = mjd2djd(mjd)
            sun.compute(Observatory)
            moon.compute(Observatory)
            pre_calced = sm.returnSunMoon(mjd)

            self.assertLess(np.abs(pre_calced['sunAlt']-sun.alt), arcmin_places)
            sun_dist = haversine(sun.ra, sun.dec, pre_calced['sunRA'], pre_calced['sunDec'])
            self.assertAlmostEqual(sun_dist, 0., places=arcmin_places)

            self.assertLess(np.abs(pre_calced['moonAlt']-moon.alt), arcmin_places)
            moon_dist = haversine(moon.ra, moon.dec, pre_calced['moonRA'], pre_calced['moonDec'])
            self.assertAlmostEqual(moon_dist, 0., places=arcmin_places)

            self.assertAlmostEqual(pre_calced['moonSunSep'], moon.phase/100.*180., places=arcmin_places)

    def testSBP(self):
        """
        Check that values are similar enough
        """

        original_model = sb.SkyModel(mags=True)
        pre_calc_model = sbp.SkyModelPre()

        hpindx = np.arange(hp.nside2npix(pre_calc_model.header['nside']))
        ra, dec = utils.hpid2RaDec(pre_calc_model.header['nside'], hpindx)

        # Run through a number of mjd values
        step = 30. / 60. / 24.  # 30 minute timestep
        nmjd = 48
        mjds = np.arange(nmjd)*step + 59573.2+0.1

        # Tolerance for difference between true and interpolated airmass
        am_tol = 0.05
        am_limit = 3.

        # Where to check the magnitudes match
        mag_am_limit = 1.5
        mag_tol = 0.2  # mags

        for mjd in mjds:
            original_model.setRaDecMjd(ra, dec, mjd, degrees=True)
            if original_model.sunAlt < np.radians(-12.):
                sky1 = original_model.returnMags()
                sky2 = pre_calc_model.returnMags(mjd)
                am1 = original_model.airmass
                am2 = pre_calc_model.returnAirmass(mjd)
                good_am = np.where((am1 >= 1.) & (am1 <= am_limit))
                diff = am1[good_am] - am2[good_am]
                # Check that the interpolated airmass is close
                assert(np.max(np.abs(diff)) < am_tol)

                for filtername in sky1.keys():
                    good = np.where((am1 < mag_am_limit) & (sky2[filtername] != hp.UNSEEN) &
                                    (np.isfinite(sky1[filtername])))
                    diff = sky1[filtername][good] - sky2[filtername][good]
                    assert(np.max(np.abs(diff)) <= mag_tol)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
