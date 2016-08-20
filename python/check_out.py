import numpy as np
import healpy as hp
import matplotlib.pylab as plt



result = np.load('generated_sky.npz')

oneD = result['dict_of_lists'][()]
maps = result['sky_brightness'][()]

result.close()

