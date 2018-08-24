#!/bin/bash
set -e

FETCH_CMD="rsync -av --progress"
DATA_URL="lsst-rsync.ncsa.illinois.edu::sim/sims_skybrightness_pre"

usage() {
  cat << EOD

  Usage: $(basename "$0") [options]

  This downloads data for the sims_skybrightness_pre package.
  Default should download the 3 and 30 day files. Use -p to
  get all the healpix data. 

  Available options:
    -h          this message
    -o          Download OpSim fields data only
    -p          Download all HealPix data

EOD
}

HEALPIX=0
OPSIM_FIELDS=0
SMALL=1
# get the options
while getopts hop c; do
    case $c in
            h) usage ; exit 0 ;;
            o) OPSIM_FIELDS=1 ; HEALPIX=0 ; SMALL=0 ;;
            p) HEALPIX=1 ; OPSIM_FIELDS=0 ; SMALL=0 ;;
            \?) usage ; exit 2 ;;
    esac
done

# Copy data down from NCSA
if [ ${HEALPIX} -eq 1 ]; then
	${FETCH_CMD} ${DATA_URL}/healpix/*.npz* healpix/
  ${FETCH_CMD} ${DATA_URL}/healpix/*.npy* healpix/
fi
if [ ${OPSIM_FIELDS} -eq 1 ]; then
	${FETCH_CMD} ${DATA_URL}/opsimFields/*.npz opsimFields/
  ${FETCH_CMD} ${DATA_URL}/opsimFields/*.npy opsimFields/
fi
if [ ${SMALL} -eq 1 ]; then
  ${FETCH_CMD} ${DATA_URL}/healpix/59853*.npz healpix/
  ${FETCH_CMD} ${DATA_URL}/healpix/59853*.npy healpix/
fi
${FETCH_CMD} ${DATA_URL}/percentile_m5_maps.npz .
