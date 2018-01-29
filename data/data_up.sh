#!/bin/bash

# Copy data up to NCSA
rsync -av --progress healpix/*.npz* lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/healpix/"
rsync -av --progress healpix/*.npy* lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/healpix/"
rsync -av --progress opsimFields/*.npz lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/opsimFields/"
rsync -av --progress opsimFields/*.npy lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/opsimFields/"
rsync -av --progress percentile_m5_maps.npz lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/"
