#! /bin/bash
source ~/.profile

# Save arguments to file
echo "faceter.sh $@" > faceter-$(date +%s).last

set -e
set -x


# Parse parameters
LONG='size:,center:,scale:,oversample:,channelsout:,threshold:,maxsize:,minsize:,datacolumn:,calopts:,wscleanopts:,weight:,noprimarybeam,mwapath:,subroutines:,facetid:'
OPTS=$(getopt --options '' --longoptions $LONG --name "$0" -- "$@")
eval set -- "$OPTS"

# Default values
size=7500
oversample=2
channelsout=6
threshold=4
maxsize=4
minsize=0.5
datacolumn=CORRECTED
calopts=""
wscleanopts=""
weight="briggs 0.5"
noprimarybeam=false
mwapath=""

while true; do
  case "$1" in
    --size )
      size="$2"
      shift 2
      ;;
    --center )
      center="$2"
      shift 2
      ;;
    --scale )
      scale="$2"
      shift 2
      ;;
    --oversample )
      oversample="$2"
      shift 2
      ;;
    --channelsout )
      channelsout="$2"
      shift 2
      ;;
    --threshold )
      threshold="$2"
      shift 2
      ;;
    --maxsize )
      maxsize="$2"
      shift 2
      ;;
    --minsize )
      minsize="$2"
      shift 2
      ;;
    --datacolumn )
      datacolumn="$2"
      shift 2
      ;;
    --calopts )
      calopts="$2"
      shift 2
      ;;
    --wscleanopts )
      wscleanopts="$2"
      shift 2
      ;;
    --weight )
      weight="$2"
      shift 2
      ;;
    --noprimarybeam )
      noprimarybeam=true
      shift 1
      ;;
    --mwapath )
      mwapath="-mwa-path $2"
      shift 2
      ;;
    --subroutines )
      IFS=','
      subroutines=( $2 )
      unset IFS
      shift 2
      ;;
    --facetid )
      facetid="$2"
      shift 2
      ;;
    -- )
      shift
      break
      ;;
    * )
      break;
      ;;
  esac
done

# Check for required parameters
if [[ -z $scale || -z $center ]]; then
  echo "Missing required parameter"
  exit 1
fi

# Check at least one mset is provided
if [[ $# < 1 ]]; then
  echo "No msets provided"
  exit 1
fi

if [[ $noprimarybeam = false ]]; then
  wscleanopts="$wscleanopts $mwapath -grid-with-beam"
  pb="-pb"
fi

# Load msets into arrays
msets=(); splits=(); facets=()
for i in $(seq 1 $#); do
  splits+=( mset-${i}.ms )
  facets+=( facet-mset-${i}.ms )
  msets+=( ${!i} )
done

echo "Processing msets: ${msets[@]}"

setup () {
  local mset
  for i in ${!msets[@]}; do
    mset=${msets[$i]}
    split=${splits[$i]}

    # Set phase center
    chgcentre $mset $center

    # Split out mset
    rm -r ${split} || true
    echo "split(vis='${mset}', outputvis='${split}', datacolumn='${datacolumn}')" | ~/casa/bin/casa -nologger --agg --nogui -c
  done
}

fullsky () {
  # First do full-sky image
  wsclean \
    -name fullsky \
    -size $size $size \
    -scale $scale \
    -mgain 0.8 \
    -niter 9999999 \
    -nmiter 12 \
    -pol i \
    -auto-threshold 3 \
    -channels-out $channelsout \
    -fit-spectral-pol 2 \
    -join-channels \
    -weight $weight \
    -parallel-deconvolution 1024 \
    -use-idg \
    -idg-mode hybrid \
    -data-column DATA \
    $wscleanopts \
    ${splits[@]}

  local mset
  for mset in ${splits[@]}; do
    # Save the data
    columnabacus.py ${mset}::DATA_SAVED = ${mset}::DATA

    # Save the model
    columnabacus.py  ${mset}::MODEL_SAVED = ${mset}::MODEL_DATA
  done
}

create_facets () {
  rm facet_centers.pkl dists.npy facet-*-model.fits || true
  local i
  for i in $(seq 0 $((channelsout - 1)) ); do
    create_facets.py --threshold $threshold --min $minsize --max $maxsize --channel 000${i} --image fullsky-MFS-image.fits --model fullsky-000${i}-model${pb}.fits
  done
}

prepare_cols() {
  local mset
  for mset in ${splits[@]}; do
    chgcentre $mset $center

    # Set corrected_data as residuals
    columnabacus.py ${mset}::CORRECTED_DATA = ${mset}::DATA_SAVED - ${mset}::MODEL_SAVED

    # Reset data
    columnabacus.py ${mset}::DATA = ${mset}::DATA_SAVED
  done
}

pre_predict() {
  # Predict facet into MODEL_DATA column
  local mset
  for mset in ${splits[@]}; do
    chgcentre $mset $center
  done

  wsclean \
    -name facet-${facetid} \
    -size $size $size  \
    -scale $scale \
    -pol i \
    -channels-out $channelsout \
    -weight $weight \
    -use-idg \
    -idg-mode hybrid \
    -predict \
    ${splits[@]}
}

split_facet() {
  for i in ${!splits[@]}; do
    local mset=${splits[$i]}
    local facet=${facets[$i]}

    # Add prediction to residuals
    columnabacus.py ${mset}::CORRECTED_DATA = ${mset}::CORRECTED_DATA + ${mset}::MODEL_DATA

    # Split out averaged ms
    rm -r ${facet} || true
    # TODO calculate width based on maximum angular size of facet
    chgcentre ${mset} $facetcenter
    echo "split(vis='${mset}', outputvis='${facet}', datacolumn='CORRECTED', width=8)" | ~/casa/bin/casa -nologger --agg --nogui -c
  done
}

set_facet_vars() {
  local model="facet-${1}-0000-model.fits"
  facetcenter=$(fitsheader -k FCTCEN -t ascii.csv $model | tail -n 1 | cut -d ',' -f 4)
  finescale=$(echo "scale=6; $scale / $oversample" | bc)
  local maxangulardist=$(fitsheader -k FCTMAX -t ascii.csv $model | tail -n 1 | cut -d ',' -f 4)
  pixels=$(echo "(($maxangulardist / $finescale) * 2.2) / 1" | bc)  # divide by 1 to round to integer
}

image_before() {
  # Image facet
  wsclean \
    -name facet-${facetid}-before \
    -size $pixels $pixels \
    -scale $finescale \
    -mgain 0.8 \
    -niter 9999999 \
    -multiscale \
    -nmiter 12 \
    -pol i \
    -auto-mask 3 \
    -auto-threshold 1 \
    -channels-out $channelsout \
    -fit-spectral-pol 2 \
    -join-channels \
    -weight $weight \
    -data-column DATA \
    -use-idg \
    -idg-mode hybrid \
    $wscleanopts \
    ${facets[@]}
}

selfcal() {
  # Selfcalibrate each input mset based on the jointly imaged model
  local i
  for i in ${!facets[@]}; do
    local facet=${facets[$i]}
    calibrate -datacolumn DATA -a 5e-6 1e-8 $calopts $facet solutions-facet-${facetid}-${i}.bin
    applysolution.py --src DATA --dest CORRECTED_DATA $facet solutions-facet-${facetid}-${i}.bin
  done

  # Create clean mask based on first imaging round
  create_mask.py --max $maxsize --facetid $facetid --image facet-${facetid}-before-MFS-image.fits
}

post_predict() {
  # Image calibrated ms and build improved model
  wsclean \
    -name facet-${facetid}-after \
    -size $pixels $pixels \
    -scale $finescale \
    -mgain 0.8 \
    -niter 9999999 \
    -multiscale \
    -nmiter 12 \
    -pol i \
    -auto-mask 3 \
    -auto-threshold 1 \
    -channels-out $channelsout \
    -fit-spectral-pol 2 \
    -join-channels \
    -weight $weight \
    -data-column CORRECTED_DATA \
    -use-idg \
    -idg-mode hybrid \
    -fits-mask facet-${facetid}-mask.fits \
    $wscleanopts \
    ${facets[@]}

  # Predict corrected facet into MODEL_DATA column
  local mset
  for mset in ${splits[@]}; do
    chgcentre $mset $facetcenter
  done
  wsclean \
    -name facet-${facetid}-after \
    -size $pixels $pixels \
    -scale $finescale \
    -pol i \
    -channels-out $channelsout \
    -weight $weight \
    -use-idg \
    -idg-mode hybrid \
    -predict \
    ${splits[@]}
}

subtract() {
  # Update residuals and data of each mset with new calibration and model data
  local i
  local mset
  for i in ${!splits[@]}; do
    mset=${splits[$i]}

    chgcentre $mset $center
    columnabacus.py ${mset}::DATA = ${mset}::DATA + ${mset}::MODEL_DATA
    applysolution.py --reverse --src MODEL_DATA --dest MODEL_DATA ${mset} solutions-facet-${facetid}-${i}.bin
    columnabacus.py ${mset}::DATA = ${mset}::DATA - ${mset}::MODEL_DATA

    # Update residuals
    columnabacus.py ${mset}::CORRECTED_DATA = ${mset}::CORRECTED_DATA - ${mset}::MODEL_DATA
  done
}

fullsky_corrected() {
  # Do fullsky image of corrected data
  wsclean \
    -name fullsky-corrected \
    -size $size $size \
    -scale $scale \
    -mgain 0.8 \
    -niter 9999999 \
    -nmiter 12 \
    -pol i \
    -auto-threshold 3 \
    -channels-out $channelsout \
    -fit-spectral-pol 2 \
    -join-channels \
    -weight $weight \
    -parallel-deconvolution 1024 \
    -use-idg \
    -idg-mode hybrid \
    -data-column DATA \
    $wscleanopts \
    ${splits[@]}
}

# Run specified subroutine if it exists
if [[ -n $subroutines ]]; then
  if [[ -n $facetid ]]; then
    set_facet_vars $facetid
  fi

  for subroutine in ${subroutines[@]}; do
    $subroutine
  done

  exit 0
fi

# Otherwise, run full faceting algorithm...
setup
fullsky
create_facets
prepare_cols

for model in facet-?-0000-model.fits facet-??-0000-model.fits; do
  if [[ ! -f $model ]]; then
    continue
  fi

  facetid=$(echo $model | cut -d '-' -f 2)
  echo "Processing facet $facetid"

  # Set global variables for this facet (facetcenter, finescale, pixels)
  set_facet_vars $facetid

  # Predict facet into model column
  pre_predict

  # Phase rotate on facet center and split into frequency averaged msets ("splits")
  split_facet

  # Image facet before (populate model column)
  image_before

  # Selfcal (and create clean mask)
  selfcal

  # Reimage after calibration and predict new model into mset
  post_predict

  # Subtract out and update new facet model
  subtract

done

fullsky_corrected
