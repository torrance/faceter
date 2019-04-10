#! /bin/bash

set -e
set -x

obsid=$1
if [[ -n $2 ]]; then
  workingdir=$2
else
  workingdir=.
fi

center=$(pointing.py ${workingdir}/${obsid}.metafits)

setup () {
  # Set phase center
  # Do this before splitting with casa in case we are set to minw
  chgcentre ${workingdir}/${obsid}.ms $center

  # Split out mset
  rm -r faceter.ms || true
  echo "split(vis='${workingdir}/${obsid}.ms', outputvis='faceter.ms', datacolumn='DATA')" | ~/casa/bin/casa -nologger --agg --nogui -c
}

fullsky () {
  # First do full-sky image
  scale=$(echo "scale=6; 0.6 / $(getchan.py ${workingdir}/${obsid}.metafits)" | bc)
  wsclean \
    -name $obsid-fullsky \
    -size 7500 7500 \
    -scale $scale \
    -mgain 0.8 \
    -niter 9999999 \
    -nmiter 12 \
    -pol i \
    -auto-threshold 3 \
    -channels-out 6 \
    -fit-spectral-pol 2 \
    -join-channels \
    -weight briggs 0.5 \
    -parallel-deconvolution 1024 \
    -use-idg \
    -idg-mode hybrid \
    -data-column DATA \
    -temp-dir /tmp \
    faceter.ms

   # Save the skymodel
   ./columnabacus.py faceter.ms::FULLSKYMODEL_DATA = faceter.ms::MODEL_DATA

   # Save the data
   ./columnabacus.py faceter.ms::DATA_SAVED = faceter.ms::DATA

   # Set corrected_data as residuals
   ./columnabacus.py faceter.ms::CORRECTED_DATA = faceter.ms::DATA - faceter.ms::MODEL_DATA
}

create_facets () {
  rm facet_centers.pkl dists.npy facet-*-model.fits || true
  local i
  for i in {0..5}; do
    ./create_facets.py --threshold 8 --min 0.5 --max 4 --channel 000${i} --image ${obsid}-fullsky-MFS-image.fits --model ${obsid}-fullsky-000${i}-model.fits
  done
}

setup
fullsky
create_facets

./columnabacus.py faceter.ms::DATA = faceter.ms::DATA_SAVED
./columnabacus.py faceter.ms::CORRECTED_DATA = faceter.ms::DATA - faceter.ms::FULLSKYMODEL_DATA

for facet in facet-?-0000-model.fits facet-??-0000-model.fits; do
  name=$(basename $facet -0000-model.fits)
  facetid=$(echo $name | cut -d '-' -f 2)
  echo "Processing $name"

  # Predict facet into MODEL_DATA column
  chgcentre faceter.ms $center
  scale=$(echo "scale=6; 0.6 / $(getchan.py ${workingdir}/${obsid}.metafits)" | bc)
  wsclean \
    -name $name \
    -size 7500 7500  \
    -scale $scale \
    -pol i \
    -channels-out 6 \
    -weight briggs 0.5 \
    -use-idg \
    -idg-mode hybrid \
    -predict \
    -temp-dir /tmp \
    faceter.ms

  # Add prediction to residuals
  ./columnabacus.py faceter.ms::CORRECTED_DATA = faceter.ms::CORRECTED_DATA + faceter.ms::MODEL_DATA

  # Split out averaged ms
  rm -r ${name}.ms || true
  # TODO calculate width based on maximum angular size of facet
  facetcenter=$(fitsheader -k FCTCEN -t ascii.csv ${name}-0000-model.fits | tail -n 1 | cut -d ',' -f 4)
  chgcentre faceter.ms $facetcenter
  echo "split(vis='faceter.ms', outputvis='${name}.ms', datacolumn='CORRECTED', width=8)" | ~/casa/bin/casa -nologger --agg --nogui -c

  # Calculate size
  scale=$(echo "scale=6; 0.3 / $(getchan.py ${workingdir}/${obsid}.metafits)" | bc)
  maxangulardist=$(fitsheader -k FCTMAX -t ascii.csv ${name}-0000-model.fits | tail -n 1 | cut -d ',' -f 4)
  pixels=$(echo "(($maxangulardist / $scale) * 2.2) / 1" | bc)  # divide by 1 to round to integer

  # Image facet
  wsclean \
    -name ${name}-before \
    -size $pixels $pixels \
    -scale $scale \
    -mgain 0.8 \
    -niter 9999999 \
    -multiscale \
    -nmiter 12 \
    -pol i \
    -auto-mask 3 \
    -auto-threshold 1 \
    -channels-out 6 \
    -fit-spectral-pol 2 \
    -join-channels \
    -weight briggs 0.5 \
    -data-column DATA \
    -use-idg \
    -idg-mode hybrid \
    -temp-dir /tmp \
    ${name}.ms

  # Self calibrate
  calibrate -datacolumn DATA -a 5e-6 1e-8 ${name}.ms solutions-${name}.bin
  ./applysolution.py --src DATA --dest CORRECTED_DATA  ${name}.ms solutions-${name}.bin

  # Create clean mask based on first imaging round
  ./create_mask.py --max 4 --facetid $facetid --image ${name}-before-MFS-image.fits

  # Image calibrated ms and build improved model
  wsclean \
    -name ${name}-after \
    -size $pixels $pixels \
    -scale $scale \
    -mgain 0.8 \
    -niter 9999999 \
    -multiscale \
    -nmiter 12 \
    -pol i \
    -auto-mask 3 \
    -auto-threshold 1 \
    -channels-out 6 \
    -fit-spectral-pol 2 \
    -join-channels \
    -weight briggs 0.5 \
    -data-column CORRECTED_DATA \
    -use-idg \
    -idg-mode hybrid \
    -fits-mask facet-${facetid}-mask.fits \
    -temp-dir /tmp \
    ${name}.ms

  # Predict corrected facet into MODEL_DATA column
  chgcentre faceter.ms $facetcenter
  wsclean \
    -name ${name}-after \
    -size $pixels $pixels \
    -scale $scale \
    -pol i \
    -channels-out 6 \
    -weight briggs 0.5 \
    -use-idg \
    -idg-mode hybrid \
    -predict \
    -temp-dir /tmp \
    faceter.ms

  # Add corrected model to fullskymodeal
  # final = correctedmodel + residuals
  #       = (fullskymodel - uncorrected_facet + corrected_facet) + residuals
  chgcentre faceter.ms $center
  ./columnabacus.py faceter.ms::DATA = faceter.ms::DATA + faceter.ms::MODEL_DATA
  ./applysolution.py --reverse --src MODEL_DATA --dest MODEL_DATA faceter.ms solutions-${name}.bin
  ./columnabacus.py faceter.ms::DATA = faceter.ms::DATA - faceter.ms::MODEL_DATA

  # Update residuals
  ./columnabacus.py faceter.ms::CORRECTED_DATA = faceter.ms::CORRECTED_DATA - faceter.ms::MODEL_DATA

done

