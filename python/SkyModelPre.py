import numpy as np
import glob

class SkyModelPre(object):

    def __init__(self, filename, data_path=None):
        """
        
        """

        # XXX.  Make things intelligent so that it takes a path and knows how to read 
        # filenames for the mjds that are inside and updates based on that

        self.info =None
        self.sb = None

        if data_path is None:
            # data_path = 
            pass
        self.files = glob.glob(os.path.join(data_path), '*.npz')

        data = np.load(filename)
        self.info = data['dict_of_lists'][()]
        self.sb = data['sky_brightness'][()]
        data.close()
        self.filter_names = self.sb.keys()

    def _load_data(self, mjd, data_path=None):
        """
        Load up the .npz file to interpolate things
        """

        # Figure out which file to load.
        filename = None

        data = np.load(filename)
        self.info = data['dict_of_lists'][()]
        self.sb = data['sky_brightness'][()]
        data.close()
        self.filter_names = self.sb.keys()

    def full_sky(self, mjd, apply_mask=False):
        """
        return a full sky map for the input mjd
        """

        if self.info is None:
            self._load_data(mjd)

        left = np.searchsorted(self.info['mjds'], mjd)-1
        right = left+1

        # Make sure the mjd is in the range
        if (left == 0) | (right >= np.size(self.info['mjds'])):
            raise ValueError("poop")

        baseline = self.info['mjds'][right] - self.info['mjds'][left]
        wterm = (mjd - self.info['mjds'][left])/baseline
        w1 = (1. - wterm)
        w2 = wterm
        sbs = {}
        for filter_name in self.filter_names:
            sbs[filter_name] = self.sb[filter_name][left, :] * w1 + self.sb[filter_name][right, :] * w2
            # XXX just do a quick nearest neighbor return
            #sbs[filter_name] = self.sb[filter_name][left, :]
        return sbs


