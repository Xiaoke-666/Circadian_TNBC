#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH --array=1-16%8
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -J cellbender

SAMPLE=$(< "../samples.txt" awk -v TASK=${SLURM_ARRAY_TASK_ID} 'NR == TASK { print $1 }')
Indir="./001_raw_cellranger_counts"
OUTSDIR="./002_Cellbender_counts"
H5AD=$Indir/"${SAMPLE}.raw_feature_bc_matrix.h5"
OUTPUT=$OUTSDIR/"${SAMPLE}.raw.cellbender.h5"
SUMMARY=${Indir}/"${SAMPLE}.summary.csv"

NEXP=$(cut -d , -f4 "${SUMMARY}" | tail -n1)
echo $NEXP
source activate cellbender
cellbender remove-background --input "${H5AD}" --output "${OUTPUT}" --low-count-threshold 5  --expected-cells "${NEXP}" 
