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

if [[ -f facet_order ]]; then
  facet_order=$(cat facet_order)
else
  facet_order="1 2 3 4 5 6 7 8 9"
fi

for facetid in $facet_order; do
  jobid=$(sbatch $gpuq_params -t 16:00:00 -d afterok:${jobid} \
          faceter.sh --subroutines pre_predict,split_facet,image_before,selfcal,post_predict,subtract \
                     --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)
done

jobid=$(sbatch $gpuq_params -d afterok:${jobid} -- \
        faceter.sh --subroutines fullsky_corrected --facetid $facetid "${params[@]}" ${msets[@]} | cut -d ' ' -f 4)

