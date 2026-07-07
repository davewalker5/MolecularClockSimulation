#!/usr/bin/env bash

if (( $# < 2 )); then
    scriptname=$(basename -- "$0")
    echo Usage: $scriptname /path/to/FASTA/file /path/to/output/folder [Hamming|Proportional]
    exit 1
fi

export PROJECT_ROOT=$( cd "$(dirname "$0")/.." ; pwd -P )
cd "$PROJECT_ROOT"

export PYTHONPATH="$PROJECT_ROOT/src"

type="${3:-hamming}"
type=$(printf '%s' "$type" | tr '[:upper:]' '[:lower:]')

python -m distancematrix --input "$1" --output "$2" --distance-type "$type"
