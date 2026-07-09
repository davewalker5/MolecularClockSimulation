#!/usr/bin/env bash

if (( $# != 1 )); then
    scriptname=$(basename -- "$0")
    echo "Usage: $scriptname /path/to/config/file"
    exit 1
fi

export PROJECT_ROOT=$( cd "$(dirname "$0")/.." ; pwd -P )
cd "$PROJECT_ROOT"

. venv/bin/activate

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
original_tree="$output_folder/true_tree.newick"

# Simulate the sequences and generate the terminal sequence file
python -m "${model}clock" --config "$1"

for distance in "${DISTANCE_TYPES[@]}"; do
    # Calculate the distance matrix using the current substitution model
    python -m distancematrix --input "$terminal_sequences" --output "$output_folder" --distance-type "$distance"

    # Reconstruct the phylogenetic tree
    python -m phylogeny --input "$output_folder/distance_matrix_${distance}.json"

    # Create the comparison image
    reconstructed_tree="$output_folder/upgma_${distance}.newick"
    comparison_image="$output_folder/upgma_${distance}.png"
    python -m treecomparison \
        --source "$original_tree" \
        --reconstructed "$reconstructed_tree" \
        --output "$comparison_image"
done
