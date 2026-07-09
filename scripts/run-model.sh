#!/usr/bin/env bash

if (( $# != 1 )); then
    scriptname=$(basename -- "$0")
    echo "Usage: $scriptname /path/to/config/file"
    exit 1
fi

export PROJECT_ROOT=$( cd "$(dirname "$0")/.." ; pwd -P )
cd "$PROJECT_ROOT"

export PYTHONPATH="$PROJECT_ROOT/src"

# Determine which model to run from the input file and validate it
model=$(sed -n 's/.*"clock_model"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$1")
case "$model" in
    strict|relaxed)
        ;;
    *)
        echo "Invalid model: $model"
        exit 1
        ;;
esac

# Define available distance calculations/substitution models
declare -a DISTANCE_TYPES=(
    hamming
    proportional
    jc69
    k80
    f81
    hky85
)

# Get the path to the output folder and the terminal sequences file
stem=$(basename "$1")
stem=${stem%.*}
output_folder="$PROJECT_ROOT/data/output/$stem"
terminal_sequences="$output_folder/terminal_sequences.fasta"

# Run the steps in the model
python -m "${model}clock" --config "$1"

for distance in "${DISTANCE_TYPES[@]}"; do
    python -m distancematrix --input "$terminal_sequences" --output "$output_folder" --distance-type "$distance"
done
