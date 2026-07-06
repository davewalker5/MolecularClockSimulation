#!/usr/bin/env bash

export PROJECT_ROOT=$( cd "$(dirname "$0")/.." ; pwd -P )
cd "$PROJECT_ROOT"

streamlit run src/strictclock/explorer.py
