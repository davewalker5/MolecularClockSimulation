[![GitHub issues](https://img.shields.io/github/issues/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/issues)
[![Releases](https://img.shields.io/github/v/release/davewalker5/MolecularClockSimulation.svg?include_prereleases)](https://github.com/davewalker5/MolecularClockSimulation/releases)
[![License](https://img.shields.io/badge/License-mit-blue.svg)](https://github.com/davewalker5/MolecularClockSimulation/blob/main/LICENSE)
[![Language](https://img.shields.io/badge/language-python-blue.svg)](https://www.python.org)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/)

# Molecular Clock Simulation

<img src="https://github.com/davewalker5/MolecularClockSimulation/blob/main/diagrams/relaxed-clock-explorer.png" width="100%">

A Python toolkit for simulating sequence evolution, exploring molecular clock models and implementing phylogenetic reconstruction algorithms from first principles.

## Overview

Molecular clocks use DNA sequence variation to estimate evolutionary relationships and divergence times between species.

This project explores the computational foundations of molecular clock analysis by building the core algorithms from first principles. It includes simulation engines for generating synthetic evolutionary datasets, analytical tools for measuring genetic distance, interactive explorers for investigating molecular clock models, and implementations of phylogenetic reconstruction algorithms.

The repository combines reusable simulation engines, command-line interfaces for generating synthetic datasets, and interactive explorers for experimenting with both strict and relaxed molecular clock models.

The project is intended both as an educational resource and as a computational laboratory for investigating phylogenetic algorithms.

## Current Features

### Strict Molecular Clock Simulation

The core simulator implements a strict molecular clock capable of generating complete synthetic evolutionary datasets.

Features include:

- Random ancestral DNA sequence generation
- Rooted phylogenetic tree generation
- Balanced and random branching models
- Ultrametric tree calibration
- DNA sequence evolution along each branch
- Complete evolutionary history tracking
- Export of simulated datasets in FASTA, Newick and JSON formats

These datasets provide known ground truth against which future phylogenetic reconstruction algorithms can be evaluated.

### Relaxed Molecular Clock Simulation

A second simulation engine implements a simple autocorrelated relaxed molecular clock.

Unlike the strict clock model, individual lineages evolve at different mutation rates while preserving a complete known evolutionary history.

Features include:

- Independent relaxed clock simulation engine
- Lineage-specific mutation rates
- Autocorrelated rate inheritance
- Time and genetic distance tracking
- Export in FASTA, Newick and JSON formats
- Complete mutation and lineage-rate history

The relaxed clock engine complements the existing strict clock simulator, providing synthetic datasets for investigating the effects of rate variation on downstream phylogenetic analysis.

### Distance Matrix Calculator

The project includes a reusable distance matrix calculator for analysing aligned DNA sequences.

The calculator computes pairwise genetic distances between every sequence in a FASTA file and produces symmetric distance matrices suitable for downstream phylogenetic reconstruction. Two distance calculations are currently supported:

- Hamming distance
- Proportional distance (p-distance)

Matrices are exported in both CSV and JSON formats, allowing them to be inspected directly or consumed programmatically by subsequent algorithms.

The calculator is implemented as a reusable Python module with a lightweight command-line interface and forms the foundation for future phylogenetic reconstruction methods, beginning with UPGMA.

### Molecular Clock Explorers

The project includes interactive Streamlit applications built directly on top of both simulation engines.

#### Strict Molecular Clock Explorer

The **Strict Molecular Clock Explorer** provides an interactive interface for experimenting with the strict molecular clock model.

Simulation parameters can be adjusted before generating a new evolutionary history and visualising the resulting ultrametric phylogenetic tree.

The explorer provides an intuitive way to investigate how sequence length, mutation rate and tree topology influence the simulated dataset.

#### Relaxed Molecular Clock Explorer

The **Relaxed Molecular Clock Explorer** extends the same interface to the relaxed molecular clock model.

In addition to the standard simulation controls, it introduces lineage-specific mutation rates and allows the effects of evolutionary rate variation to be explored interactively.

The explorer visualises branch-specific genetic change, observed substitutions and simulation summary statistics, providing an intuitive way to investigate how relaxed molecular clocks differ from the strict clock model.

Both explorers are built on the underlying command-line simulation engines, which remain fully supported for reproducible dataset generation and downstream analysis.

## Getting Started

Further documentation, including the documentation on intended field usage and how to run the application, is available in the project [Wiki](https://github.com/davewalker5/MolecularClockSimulation/wiki).

## Feedback

To file issues or suggestions, please use the [Issues](https://github.com/davewalker5/MolecularClockSimulation/issues) page for this project on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
