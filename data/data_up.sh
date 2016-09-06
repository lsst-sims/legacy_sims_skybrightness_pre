#!/bin/bash

# Copy data up to NCSA
rsync -av --progress healpix/*.npz lsst-dev.ncsa.illinois.edu:"/lsst/sim-data/sims_skybrightness_pre/healpix/"
rsync -av --progress opsimFields/*.npz lsst-dev.ncsa.illinois.edu:"/lsst/sim-data/sims_skybrightness_pre/opsimFields/"
