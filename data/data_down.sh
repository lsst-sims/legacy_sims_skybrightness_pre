#!/bin/bash

# Copy data down from NCSA
rsync -av --progress lsst-rsync.ncsa.illinois.edu::sim/sims_skybrightness_pre/healpix/*.npz healpix/
rsync -av --progress lsst-rsync.ncsa.illinois.edu::sim/sims_skybrightness_pre/opsimFields/*.npz opsimFields/
