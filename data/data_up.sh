#!/bin/bash

# Copy data up to NCSA
rsync -avz --progress healpix/*.npz* lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/healpix/"
rsync -avz --progress healpix/*.npy* lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/healpix/"
rsync -avz --progress opsimFields/*.npz lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/opsimFields/"
rsync -avz --progress opsimFields/*.npy lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/opsimFields/"
rsync -avz --progress percentile_m5_maps.npz lsst-dev01.ncsa.illinois.edu:"/datasets/public_html/sim-data/sims_skybrightness_pre/"
