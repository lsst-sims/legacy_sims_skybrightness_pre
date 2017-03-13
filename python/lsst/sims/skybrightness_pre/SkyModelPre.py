import numpy as np
import glob
import os
import healpy as hp
from lsst.utils import getPackageDir
import warnings
from lsst.sims.utils import haversine

__all__ = ['SkyModelPre']


class SkyModelPre(object):
    """
    Load pre-computed sky brighntess maps for the LSST site and use them to interpolate to
    arbitrary dates.
    """

    def __init__(self, data_path=None, opsimFields=False, preload=True, speedLoad=False):

        self.info = None
        self.sb = None
        self.opsimFields = opsimFields

        # Look in default location for .npz files to load
        if 'SIMS_SKYBRIGHTNESS_DATA' in os.environ:
            data_dir = os.environ['SIMS_SKYBRIGHTNESS_DATA']
        else:
            data_dir = os.path.join(getPackageDir('sims_skybrightness_pre'), 'data')

        if data_path is None:
            if opsimFields:
                data_path = os.path.join(data_dir, 'opsimFields')
            else:
                data_path = os.path.join(data_dir, 'healpix')

        self.files = glob.glob(os.path.join(data_path, '*.npz*'))
        if len(self.files) == 0:
            errmssg = 'Failed to find pre-computed .npz files. '
            errmssg += 'Copy data from NCSA with sims_skybrightness_pre/data/data_down.sh \n'
            errmssg += 'or build by running sims_skybrightness_pre/data/generate_sky.py'
            raise ValueError(errmssg)
        mjd_left = []
        mjd_right = []
        # Expect filenames of the form mjd1_mjd2.npz, e.g., 59632.155_59633.2.npz
        for filename in self.files:
            temp = os.path.split(filename)[-1].replace('.npz', '').split('_')
            mjd_left.append(float(temp[0]))
            mjd_right.append(float(temp[1]))

        self.mjd_left = np.array(mjd_left)
        self.mjd_right = np.array(mjd_right)

        # Go ahead and load the first one by default
        if speedLoad:
            self._load_data(59580., filename=os.path.join(data_dir, 'healpix/small_example.npz_small'))
        else:
            if preload:
                self._load_data(self.mjd_left[0])
            else:
                self.loaded_range = np.array([-1])

    def _load_data(self, mjd, filename=None):
        """
        Load up the .npz file to interpolate things
        """

        # Figure out which file to load.
        file_indx = np.where((mjd >= self.mjd_left) & (mjd <= self.mjd_right))[0]
        if np.size(file_indx) == 0:
            raise ValueError('MJD = %f is out of range for the files found (%f-%f)' % (mjd,
                                                                                       self.mjd_left.min(),
                                                                                       self.mjd_right.max()))
        if filename is None:
            filename = self.files[file_indx.min()]

        data = np.load(filename)
        self.info = data['dict_of_lists'][()]
        self.sb = data['sky_brightness'][()]
        self.header = data['header'][()]
        data.close()
        self.filter_names = self.sb.keys()

        if not self.opsimFields:
            self.nside = hp.npix2nside(self.sb[self.filter_names[0]][0, :].size)
        self.loaded_range = np.array([self.mjd_left[file_indx], self.mjd_right[file_indx]])
        if filename is not None:
            self.loaded_range = np.array([self.info['mjds'].min(), self.info['mjds'].max()])

    def returnSunMoon(self, mjd):
        """
        Return dictionary with the interpolated positions for sun and moon

        Parameters
        ----------
        mjd : float
           Modified Julian Date to interpolate to

        Returns
        -------
        sunMoon : dict
            Dict with keys for the sun and moon RA and Dec and the
            mooon-sun separation.
        """

        keys = ['sunAlts', 'moonAlts', 'moonRAs', 'moonDecs', 'sunRAs',
                'sunDecs', 'moonSunSep']

        if (mjd < self.loaded_range.min() or (mjd > self.loaded_range.max())):
            self._load_data(mjd)

        left = np.searchsorted(self.info['mjds'], mjd)-1
        right = left+1

        # If we are out of bounds
        if right >= self.info['mjds'].size:
            right -= 1
            baseline = 1.
        elif left < 0:
            left += 1
            baseline = 1.
        else:
            baseline = self.info['mjds'][right] - self.info['mjds'][left]

        wterm = (mjd - self.info['mjds'][left])/baseline
        w1 = (1. - wterm)
        w2 = wterm

        result = {}
        for key in keys:
            if key[-1] == 's':
                newkey = key[:-1]
            else:
                newkey = key
            result[newkey] = self.info[key][left] * w1 + self.info[key][right] * w2
        return result

    def returnAirmass(self, mjd, maxAM=10., indx=None, badval=hp.UNSEEN):
        """

        Parameters
        ----------
        mjd : float
            Modified Julian Date to interpolate to
        indx : List of int(s) (None)
            indices to interpolate the sky values at. Returns full sky if None. If the class was
            instatiated with opsimFields, indx is the field ID, otherwise it is the healpix ID.
        maxAM : float (10)
            The maximum airmass to return, everything above this airmass will be set to badval

        Returns
        -------
        airmass : np.array
            Array of airmass values. If the MJD is between sunrise and sunset, all values are masked.
        """
        if (mjd < self.loaded_range.min() or (mjd > self.loaded_range.max())):
            self._load_data(mjd)

        left = np.searchsorted(self.info['mjds'], mjd)-1
        right = left+1

        # If we are out of bounds
        if right >= self.info['mjds'].size:
            right -= 1
            baseline = 1.
        elif left < 0:
            left += 1
            baseline = 1.
        else:
            baseline = self.info['mjds'][right] - self.info['mjds'][left]

        if indx is None:
            result_size = self.sb[self.sb.keys()[0]][left, :].size
            indx = np.arange(result_size)
        else:
            result_size = len(indx)
        # Check if we are between sunrise/set
        if baseline > self.header['timestep_max']:
            warnings.warn('Requested MJD between sunrise and sunset, returning closest maps')
            diff = np.abs(self.info['mjds'][left:right+1]-mjd)
            closest_indx = np.array([left, right])[np.where(diff == np.min(diff))]
            airmass = self.info['airmass'][closest_indx, indx]
            mask = np.where((self.info['airmass'][closest_indx, indx].ravel() < 1.) |
                            (self.info['airmass'][closest_indx, indx].ravel() > maxAM))
            airmass = airmass.ravel()

        else:
            wterm = (mjd - self.info['mjds'][left])/baseline
            w1 = (1. - wterm)
            w2 = wterm
            airmass = self.info['airmass'][left, indx] * w1 + self.info['airmass'][right, indx] * w2
            mask = np.where((self.info['airmass'][left, indx] < 1.) |
                            (self.info['airmass'][left, indx] > maxAM) |
                            (self.info['airmass'][right, indx] < 1.) |
                            (self.info['airmass'][right, indx] > maxAM))
        airmass[mask] = badval

        return airmass

    def returnMags(self, mjd, indx=None, airmass_mask=True, planet_mask=True,
                   moon_mask=True, zenith_mask=True, badval=hp.UNSEEN,
                   filters=['u', 'g', 'r', 'i', 'z', 'y'], extrapolate=False):
        """
        Return a full sky map or individual pixels for the input mjd

        Parameters
        ----------
        mjd : float
            Modified Julian Date to interpolate to
        indx : List of int(s) (None)
            indices to interpolate the sky values at. Returns full sky if None. If the class was
            instatiated with opsimFields, indx is the field ID, otherwise it is the healpix ID.
        airmass_mask : bool (True)
            Set high (>2.5) airmass pixels to badval.
        planet_mask : bool (True)
            Set sky maps to badval near (2 degrees) bright planets.
        moon_mask : bool (True)
            Set sky maps near (10 degrees) the moon to badval.
        zenith_mask : bool (True)
            Set sky maps at high altitude (>86.5) to badval.
        badval : float (-1.6375e30)
            Mask value. Defaults to the healpy mask value.
        filters : list
            List of strings for the filters that should be returned.
        extrapolate : bool (False)
            In indx is set, extrapolate any masked pixels to be the same as the nearest non-masked
            value from the full sky map.

        Returns
        -------
        sbs : dict
            A dictionary with filter names as keys and np.arrays as values which
            hold the sky brightness maps in mag/sq arcsec.
        """
        if (mjd < self.loaded_range.min() or (mjd > self.loaded_range.max())):
            self._load_data(mjd)

        mask_rules = {'airmass': airmass_mask, 'planet': planet_mask,
                      'moon': moon_mask, 'zenith': zenith_mask}

        left = np.searchsorted(self.info['mjds'], mjd)-1
        right = left+1

        # Do full sky by default
        if indx is None:
            indx = np.arange(self.sb['r'].shape[1])
            full_sky = True
        else:
            full_sky = False

        # If we are out of bounds
        if right >= self.info['mjds'].size:
            right -= 1
            baseline = 1.
        elif left < 0:
            left += 1
            baseline = 1.
        else:
            baseline = self.info['mjds'][right] - self.info['mjds'][left]

        # Check if we are between sunrise/set
        if baseline > self.header['timestep_max']:
            warnings.warn('Requested MJD between sunrise and sunset, returning closest maps')
            diff = np.abs(self.info['mjds'][left:right+1]-mjd)
            closest_indx = np.array([left, right])[np.where(diff == np.min(diff))].min()
            sbs = {}
            for filter_name in filters:
                sbs[filter_name] = self.sb[filter_name][closest_indx, indx]
                for mask_name in mask_rules:
                    if mask_rules[mask_name]:
                        toMask = np.where(self.info[mask_name+'_masks'][closest_indx, indx])
                        sbs[filter_name][toMask] = badval
                sbs[filter_name][np.isinf(sbs[filter_name])] = badval
                sbs[filter_name][np.where(sbs[filter_name] == hp.UNSEEN)] = badval
        else:
            wterm = (mjd - self.info['mjds'][left])/baseline
            w1 = (1. - wterm)
            w2 = wterm
            sbs = {}
            for filter_name in filters:
                sbs[filter_name] = self.sb[filter_name][left, indx] * w1 + \
                    self.sb[filter_name][right, indx] * w2
                for mask_name in mask_rules:
                    if mask_rules[mask_name]:
                        toMask = np.where(self.info[mask_name+'_masks'][left, indx] |
                                          self.info[mask_name+'_masks'][right, indx] |
                                          np.isinf(sbs[filter_name]))
                        sbs[filter_name][toMask] = badval
                sbs[filter_name][np.where(sbs[filter_name] == hp.UNSEEN)] = badval
                sbs[filter_name][np.where(sbs[filter_name] == hp.UNSEEN)] = badval

        # If requested a certain pixel(s), and want to extrapolate.
        if (not full_sky) & extrapolate:
            masked_pix = False
            for filter_name in filters:
                if (badval in sbs[filter_name]) | (True in np.isnan(sbs[filter_name])):
                    masked_pix = True
            if masked_pix:
                # We have pixels that are masked that we want reasonable values for
                full_sky_sb = self.returnMags(mjd, airmass_mask=False, planet_mask=False, moon_mask=False,
                                              zenith_mask=False, filters=filters)
                good = np.where((full_sky_sb[filters[0]] != badval) & ~np.isnan(full_sky_sb[filters[0]]))[0]
                ra_full = np.radians(self.header['ra'][good])
                dec_full = np.radians(self.header['dec'][good])
                for filtername in filters:
                    full_sky_sb[filtername] = full_sky_sb[filtername][good]
                # Going to assume the masked pixels are the same in all filters
                masked_indx = np.where((sbs[filters[0]].ravel() == badval) |
                                       np.isnan(sbs[filters[0]].ravel()))[0]
                for i, mi in enumerate(masked_indx):
                    # Note, this is going to be really slow for many pixels, should use a kdtree
                    dist = haversine(np.radians(self.header['ra'][indx][i]),
                                     np.radians(self.header['dec'][indx][i]),
                                     ra_full, dec_full)
                    closest = np.where(dist == dist.min())[0]
                    for filtername in filters:
                        sbs[filtername].ravel()[mi] = np.min(full_sky_sb[filtername][closest])

        return sbs

