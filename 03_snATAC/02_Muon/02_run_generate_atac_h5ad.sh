#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p <partition>
#SBATCH -J generate_h5ad
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

# activate environment
source ~/.bashrc
conda activate <env_name>

# configurable variables
SCRIPT="<path_to_script>"
INPUT_DIR="<input_directory>"
OUTPUT_DIR="<output_directory>"
PEAK_FILE="<peak_file>"

SAMPLE="<sample_id>"

mkdir -p logs
mkdir -p "${OUTPUT_DIR}"

echo "Running sample: ${SAMPLE}"

python "${SCRIPT}" \
    -i "${INPUT_DIR}" \
    -o "${OUTPUT_DIR}" \
    -p "${PEAK_FILE}" \
    -s "${SAMPLE}"

echo "Done: ${SAMPLE}"