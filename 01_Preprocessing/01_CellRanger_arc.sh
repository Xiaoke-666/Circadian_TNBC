#!/usr/bin/env bash
#
# Run Cell Ranger ARC count using a SLURM array job
#

#SBATCH -p <partition>
#SBATCH -J cellranger_arc
#SBATCH --array=1-N%M
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail

# -----------------------------
# Configuration (EDIT THESE)
# -----------------------------
CELLRANGER_ARC="/path/to/cellranger-arc"
WORKDIR="/path/to/workdir"
SAMPLE_LIST="${WORKDIR}/sample_list.txt"

PEAKS="/path/to/peak.bed"
REFERENCE="/path/to/reference"

LOCALMEM=160

# -----------------------------
# Setup
# -----------------------------
mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

# -----------------------------
# Get sample for this task
# -----------------------------
TASK_ITEM="$(awk -v task="${SLURM_ARRAY_TASK_ID}" 'NR == task {print $1}' "${SAMPLE_LIST}")"

if [[ -z "${TASK_ITEM}" ]]; then
    echo "ERROR: no sample found for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

LIBRARIES_CSV="${WORKDIR}/${TASK_ITEM}.csv"

if [[ ! -f "${LIBRARIES_CSV}" ]]; then
    echo "ERROR: missing CSV for sample ${TASK_ITEM}" >&2
    exit 1
fi

echo "========================================"
echo "Running Cell Ranger ARC"
echo "Sample: ${TASK_ITEM}"
echo "========================================"

# -----------------------------
# Run
# -----------------------------
/usr/bin/time -v "${CELLRANGER_ARC}" count \
    --id="${TASK_ITEM}" \
    --reference="${REFERENCE}" \
    --libraries="${LIBRARIES_CSV}" \
    --peaks="${PEAKS}" \
    --localmem="${LOCALMEM}"

echo "Completed: ${TASK_ITEM}"