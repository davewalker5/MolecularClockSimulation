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

Version 0.1.0 implements a strict molecular clock simulator capable of generating complete synthetic evolutionary datasets.

The simulator:

- Generates random ancestral DNA sequences
- Constructs rooted phylogenetic tree topologies
- Calibrates ultrametric trees under a strict molecular clock
- Simulates sequence evolution along each branch
- Records the complete evolutionary history
- Exports terminal sequences and ground truth in FASTA, Newick and JSON formats

These datasets provide a controlled environment for developing and evaluating molecular clock inference algorithms.

Future development will focus on reconstructing phylogenetic trees and estimating divergence times from the simulated sequence data.

## Feedback

To file issues or suggestions, please use the [Issues](https://github.com/davewalker5/MolecularClockSimulation/issues) page for this project on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
