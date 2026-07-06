[![GitHub issues](https://img.shields.io/github/issues/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/issues)
[![Releases](https://img.shields.io/github/v/release/davewalker5/MolecularClockSimulation.svg?include_prereleases)](https://github.com/davewalker5/MolecularClockSimulation/releases)
[![License](https://img.shields.io/badge/License-mit-blue.svg)](https://github.com/davewalker5/MolecularClockSimulation/blob/main/LICENSE)
[![Language](https://img.shields.io/badge/language-python-blue.svg)](https://www.python.org)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/)

# Molecular Clock Simulation

<img src="https://github.com/davewalker5/MolecularClockSimulation/blob/main/diagrams/strict-clock-explorer.png" width="100%">

A Python toolkit for simulating sequence evolution, exploring molecular clock models and implementing phylogenetic reconstruction algorithms from first principles.

## Overview

Molecular clocks use DNA sequence variation to estimate evolutionary relationships and divergence times between species.

This project explores the computational foundations of molecular clock analysis by building the core algorithms from first principles. Rather than relying on existing phylogenetics libraries, the aim is to understand how evolutionary histories are simulated, represented, reconstructed and ultimately interpreted.

The repository combines a reusable simulation library, a command-line interface for generating synthetic datasets, and an interactive explorer for experimenting with molecular clock models.

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

### Strict Molecular Clock Explorer

Version 0.2.0 introduces the **Strict Molecular Clock Explorer**, an interactive Streamlit application built directly on top of the simulation engine.

The explorer allows simulation parameters to be adjusted interactively before generating a new evolutionary history and visualising the resulting phylogenetic tree.

The command-line interface remains fully supported and continues to provide reproducible dataset generation for downstream analysis.

## Getting Started

Further documentation, including the documentation on intended field usage and how to run the application, is available in the project [Wiki](https://github.com/davewalker5/MolecularClockSimulation/wiki).

## Feedback

To file issues or suggestions, please use the [Issues](https://github.com/davewalker5/MolecularClockSimulation/issues) page for this project on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
