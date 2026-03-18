#!/bin/bash

set -euo pipefail

JULIA_CMD="${JULIA_CMD:-julia}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RUN_SCRIPT="$REPO_ROOT/src/julia/run_phase_inference.jl"
DATA_ROOT="$REPO_ROOT/data"
RESULTS_ROOT="$REPO_ROOT/results"

DATA_INPUT="${1:-}"
OUTPUT_DIR="${2:-}"

if [ -z "$DATA_INPUT" ]; then
    echo "Usage: $0 <dataset_name_or_path> [output_dir]"
    exit 1
fi

# Auto-detect if input is just a name or full path
if [ -d "$DATA_INPUT" ]; then
    DATA_PATH="$DATA_INPUT"
elif [ -d "$DATA_ROOT/$DATA_INPUT" ]; then
    DATA_PATH="$DATA_ROOT/$DATA_INPUT"
else
    echo "Error: Cannot find dataset '$DATA_INPUT'"
    exit 1
fi

# Check if this is a direct dataset or has subdatasets
if [ -f "$DATA_PATH/expression.csv" ]; then
    # Direct dataset - run once
    DATASET_NAME=$(basename "$DATA_PATH")
    
    if [ -z "$OUTPUT_DIR" ]; then
        TIMESTAMP=$(date +%Y%m%d_%H%M%S)
        OUTPUT_DIR="$RESULTS_ROOT/${DATASET_NAME}_${TIMESTAMP}"
    fi
    
    mkdir -p "$OUTPUT_DIR"
    
    echo "Running phase inference: $DATASET_NAME"
    echo "Output: $OUTPUT_DIR"

    "$JULIA_CMD" "$RUN_SCRIPT" "$DATA_PATH" "$OUTPUT_DIR"
    echo "Completed: $OUTPUT_DIR"
else
    echo "No expression.csv found. Processing subdatasets in $DATA_PATH"
    
    PARENT_NAME=$(basename "$DATA_PATH")
    
    for SUBDIR in "$DATA_PATH"/*/; do
        [ ! -d "$SUBDIR" ] && continue
        
        SUBNAME=$(basename "$SUBDIR")
        
        if [ ! -f "$SUBDIR/expression.csv" ]; then
            echo "Skipping $SUBNAME (missing expression.csv)"
            continue
        fi
        
        if [ -z "$OUTPUT_DIR" ]; then
            TIMESTAMP=$(date +%Y%m%d_%H%M%S)
            SUB_OUTPUT="$RESULTS_ROOT/${PARENT_NAME}_${SUBNAME}_${TIMESTAMP}"
        else
            SUB_OUTPUT="${OUTPUT_DIR}/${SUBNAME}"
        fi
        
        mkdir -p "$SUB_OUTPUT"
        
        echo "Running phase inference: $PARENT_NAME/$SUBNAME"
        echo "Output: $SUB_OUTPUT"

        "$JULIA_CMD" "$RUN_SCRIPT" "${SUBDIR%/}" "$SUB_OUTPUT" && echo "Completed: $SUBNAME" || echo "Failed: $SUBNAME"
    done
    
    if [ -n "$OUTPUT_DIR" ]; then
        echo "All subdatasets processed under: $OUTPUT_DIR"
    else
        echo "All subdatasets processed under: $RESULTS_ROOT"
    fi
fi