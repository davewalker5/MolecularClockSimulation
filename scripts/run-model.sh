#!/usr/bin/env bash

if (( $# != 3 )); then
    scriptname=$(basename -- "$0")
    echo "Usage: $scriptname /path/to/config/file CLOCK_MODEL DISTANCE_MODEL"
    exit 1
fi

export PROJECT_ROOT=$( cd "$(dirname "$0")/.." ; pwd -P )
cd "$PROJECT_ROOT"

export PYTHONPATH="$PROJECT_ROOT/src"

# Determine which model to run and validate the choice
model=$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')
case "$model" in
    strict)
        module="strictclock"
        ;;
    relaxed)
        module="relaxedclock"
        ;;
    *)
        echo "Invalid model: $model"
        exit 1
        ;;
esac

# Determine which distance matrix calculation to use and validate the choice
distance=$(printf '%s' "$3" | tr '[:upper:]' '[:lower:]')
case "$distance" in
    hamming|proportional|jc69|k80|f81|hky85)
        ;;
    *)
        echo "Invalid distance calculation: $distance"
        exit 1
        ;;
esac

# Get the path to the output folder and the terminal sequences file
stem=$(basename "$1")
stem=${stem%.*}
output_folder="$PROJECT_ROOT/data/output/$stem"
terminal_sequences="$output_folder/terminal_sequences.fasta"

# Run the steps in the model
python -m "$module" --config "$1"
python -m distancematrix --input "$terminal_sequences" --output "$output_folder" --distance-type "$distance"
