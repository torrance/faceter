#! /bin/bash

params=(--center "-00h52m41.25s -24d56m10.62s" --size 1000 --scale 0.005 --oversample 1.2 --channelsout 3 --maxsize 3 --minsize 0.5 --datacolumn DATA --mwapath $BASEDIR)
workq_params="-M zeus --partition workq --account pawsey0293 -t 24:00:00 --nodes 1 --cpus_per_task 14 --mem 60G --export=BASEDIR"
gpuq_params="-M zeus --partition gpuq --partition p100 --account pawsey0293 -t 24:00:00 --nodes 1 --cpus_per_task 14 --mem 110G --export=BASEDIR"

msets=()
for i in $(seq 1 $#); do
  msets+=( ${!i} )
done

set -e
set -x

jobid=$(sbatch $workq_params -- \
        faceter.sh --subroutines prepare_cols "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)

for model in facet-?-0000-model.fits facet-??-0000-model.fits; do
  if [[ ! -f $model ]]; then
    continue
  fi

  facetid=$(echo $model | cut -d '-' -f 2)

  jobid=$(sbatch $gpu_params -d afterok:${jobid} -- \
          faceter.sh --subroutines pre_predict --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)
  jobid=$(sbatch $workq_params -d afterok:${jobid} -- \
          faceter.sh --subroutines split_facet --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)
  jobid=$(sbatch $gpuq_params-d afterok:${jobid} -- \
          faceter.sh --subroutines image_before --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)
  jobid=$(sbatch $workq_params -d afterok:${jobid} -- \
          faceter.sh --subroutines selfcal --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)
  jobid=$(sbatch $gpuq_params -d afterok:${jobid} -- \
          faceter.sh --subroutines post_predict --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)
  jobid=$(sbatch $workq_params -d afterok:${jobid} -- \
          faceter.sh --subroutines subtract --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)
done

jobid=$(sbatch $gpuq_params -d afterok:${jobid} -- \
        faceter.sh --subroutines fullsky_corrected --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)

