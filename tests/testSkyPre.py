import unittest
import lsst.sims.skybrightness as sb
import lsst.sims.skybrightness_pre as sbp
import lsst.utils.tests
import healpy as hp
import numpy as np
import lsst.sims.utils as utils


class TestSkyPre(unittest.TestCase):

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
