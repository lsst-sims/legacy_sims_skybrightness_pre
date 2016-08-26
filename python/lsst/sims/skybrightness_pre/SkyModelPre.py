import numpy as np
import glob
import os
import healpy as hp
from lsst.utils import getPackageDir

__all__ = ['SkyModelPre']


class SkyModelPre(object):
    """
    Load pre-computed sky brighntess maps for the LSST site and use them to interpolate to
    arbitrary dates.
    """

    def __init__(self, data_path=None, opsimFields=False, preload=True):

        self.info = None
        self.sb = None
        self.opsimFields = opsimFields

        # Look in default location for .npz files to load
        if data_path is None:
            if opsimFields:
                data_path = os.path.join(getPackageDir('sims_skybrightness_pre'), 'data/opsimFields')
            else:
                data_path = os.path.join(getPackageDir('sims_skybrightness_pre'), 'data/healpix')

        self.files = glob.glob(os.path.join(data_path, '*.npz'))
        mjd_left = []
        mjd_right = []
        # Expect filenames of the form mjd1_mjd2.npz, e.g., 59632.155_59633.2.npz
        for filename in self.files:
            temp = filename.replace('.npz', '').split('_')
            mjd_left.append(temp[0])
            mjd_right.append(temp[1])

        # Go ahead and load the first one by default
        if preload:
            self.mjd_left = np.array(mjd_left)
            self.mjd_right = np.array(mjd_right)
            self._load_data(mjd_left[0])

    def _load_data(self, mjd):
        """
        Load up the .npz file to interpolate things
        """

        # Figure out which file to load.
        file_indx = np.where((mjd >= self.mjd_left) & (mjd <= self.mjd_right))[0]
        if np.size(file_indx) == 0:
            raise ValueError('MJD = %f is out of range for the files found' % mjd)
        filename = self.files[file_indx.min()]

        data = np.load(filename)
        self.info = data['dict_of_lists'][()]
        self.sb = data['sky_brightness'][()]
        data.close()
        self.filter_names = self.sb.keys()

        if not self.opsimFields:
            self.nside = hp.npix2nside(self.sb[self.filter_names[0]][0, :].size)

    def returnMags(self, mjd, indx=None, apply_mask=True, badval=hp.UNSEEN,
                   filters=['u', 'g', 'r', 'i', 'z', 'y']):
        """
        Return a full sky map or individual pixels for the input mjd

        Parameters
        ----------
        mjd : float
            Modified Julian Date to interpolate to
        indx : int(s) (None)
            indices to interpolate the sky values at. Returns full sky if None. If the class was
            instatiated with opsimFields, indx is the field ID, otherwise it is the healpix ID.
        apply_mask : bool (True)
            Set sky maps to badval for regions that should be avoided (high airmass, near moon, near planets)
        badval : float (-1.6375e30)
            Mask value. Defaults to the healpy mask value.
        filters : list
            List of strings for the filters that should be returned.

        Returns
        -------
        sbs : dict
            A dictionary with filter names as keys and np.arrays as values which
            hold the sky brightness maps in mag/sq arcsec.
        """

        if (mjd < self.info['mjds'].min()) | (mjd > self.info['mjds'].max()):
            self._load_data(mjd)

        left = np.searchsorted(self.info['mjds'], mjd)-1
        right = left+1

        baseline = self.info['mjds'][right] - self.info['mjds'][left]
        wterm = (mjd - self.info['mjds'][left])/baseline
        w1 = (1. - wterm)
        w2 = wterm
        sbs = {}
        for filter_name in filters:
            if indx is None:
                sbs[filter_name] = self.sb[filter_name][left, :] * w1 + self.sb[filter_name][right, :] * w2
                if apply_mask:
                    toMask = np.where(self.info['masks'][left, :] | self.info['masks'][right, :])
                    sbs[filter_name][toMask] = badval
            else:
                sbs[filter_name] = self.sb[filter_name][left, indx] * w1 + self.sb[filter_name][right, indx] * w2
                if apply_mask:
                    toMask = np.where(self.info['masks'][left, indx] | self.info['masks'][right, indx])
                    sbs[filter_name][toMask] = badval
        return sbs

