import numpy as np
import lsst.sims.skybrightness_pre as sb

if __name__ == '__main__':
    sm = sb.SkyModelPre(preload=False, verbose=True)

    mjd0 = 59031.
    stepsize = 10. # 60
    mjds = np.arange(mjd0, mjd0+365.25*10+stepsize, stepsize)
    for mjd in mjds:
        mags = sm.returnMags(mjd)

    

"""
I think most of these are going to be loading off my spinny disk:

Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59194_59407.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59194_59407.npy
59194_59407.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59377_59590.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59377_59590.npy
59377_59590.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59560_59772.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59560_59772.npy
59560_59772.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59743_59955.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59743_59955.npy
59743_59955.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59926_60138.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/59926_60138.npy
59926_60138.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60109_60321.npz

also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60109_60321.npy
60109_60321.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60292_60504.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60292_60504.npy
60292_60504.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60475_60687.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60475_60687.npy
60475_60687.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60658_60870.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60658_60870.npy
60658_60870.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60841_61053.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/60841_61053.npy
60841_61053.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61024_61236.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61024_61236.npy
61024_61236.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61207_61419.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61207_61419.npy
61207_61419.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61390_61602.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61390_61602.npy
61390_61602.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61573_61785.npz

also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61573_61785.npy
61573_61785.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61756_61968.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61756_61968.npy
61756_61968.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61939_62151.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/61939_62151.npy
61939_62151.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/62122_62334.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/62122_62334.npy
62122_62334.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/62305_62517.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/62305_62517.npy
62305_62517.npz loaded
Loading file /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/62488_62700.npz
also loading /Users/yoachim/gitRepos/sims_skybrightness_pre/data/healpix/62488_62700.npy
62488_62700.npz loaded

real    59m15.297s
user    2m52.221s
sys     2m9.841s

# So, that took almost exactly an hour to execute. Let's do it every 10 days just to 
# be sure it's loading dominated

real    56m32.622s
user    2m48.300s
sys     2m10.474s

"""