[![GitHub issues](https://img.shields.io/github/issues/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/issues)
[![Releases](https://img.shields.io/github/v/release/davewalker5/MolecularClockSimulation.svg?include_prereleases)](https://github.com/davewalker5/MolecularClockSimulation/releases)
[![License](https://img.shields.io/badge/License-mit-blue.svg)](https://github.com/davewalker5/MolecularClockSimulation/blob/main/LICENSE)
[![Language](https://img.shields.io/badge/language-python-blue.svg)](https://www.python.org)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/)

# Molecular Clock Simulation

A Python toolkit for simulating sequence evolution and exploring molecular clock algorithms.

## Overview

Molecular clocks use DNA sequence variation to estimate evolutionary relationships and divergence times between species. This repository explores the computational techniques behind molecular clock analysis by building them from first principles.

Rather than relying on existing phylogenetics libraries, the project implements the core concepts step by step, providing an opportunity to understand how evolutionary histories can be simulated, reconstructed and analysed.

The project is intended as both a computational exploration and an educational reference.

## Current Status

Development is currently focused on the simulation framework.

The initial implementation generates synthetic datasets by:

- generating a random ancestral DNA sequence;
- constructing rooted phylogenetic trees;
- calibrating ultrametric trees under a strict molecular clock;
- simulating sequence evolution along each branch; and
- exporting datasets in standard formats for later analysis.

These simulated datasets provide a known evolutionary history ("ground truth") against which molecular clock algorithms can be evaluated.

## Planned Features

The project is expected to grow to include:

- strict molecular clock simulation;
- relaxed molecular clock models;
- alternative nucleotide substitution models;
- phylogenetic tree reconstruction;
- molecular clock calibration algorithms;
- comparative evaluation of clock methods; and
- visualisation of evolutionary histories and sequence evolution.

As with the simulation framework, future components will prioritise clarity and understanding over algorithmic complexity.

## Repository Structure

```
simulation/         Sequence evolution simulator
docs/               Documentation and design notes
examples/           Example configuration files
```

Additional components will be added as the project develops.

## Feedback

To file issues or suggestions, please use the [Issues](https://github.com/davewalker5/MolecularClockSimulation/issues) page for this project on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
