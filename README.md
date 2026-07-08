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

This project explores the computational foundations of molecular clock analysis by implementing the core algorithms from first principles. It combines reusable simulation engines, analytical tools and interactive visualisations to provide an educational framework for understanding sequence evolution and phylogenetic reconstruction.

Detailed documentation, implementation notes and usage guides are available in the project Wiki.

## Current Features

| Component                          | Description                                                                                                                                             |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Strict Molecular Clock Simulation  | Simulates sequence evolution under a strict molecular clock, producing synthetic FASTA, Newick and JSON datasets                                        |
| Relaxed Molecular Clock Simulation | Simulates lineage-specific mutation rates using a simple autocorrelated relaxed clock model                                                             |
| Distance Matrix Calculator | Generates pairwise evolutionary distance matrices from aligned DNA sequences using a selection of nucleotide substitution models |
| Strict Molecular Clock Explorer    | Interactive Streamlit application for exploring strict molecular clock simulations                                                                      |
| Relaxed Molecular Clock Explorer   | Interactive Streamlit application for investigating relaxed molecular clock behaviour                                                                   |

## Documentation

Further documentation, including the documentation on intended field usage and how to run the application, is available in the project [Wiki](https://github.com/davewalker5/MolecularClockSimulation/wiki).

## Feedback

To file issues or suggestions, please use the [Issues](https://github.com/davewalker5/MolecularClockSimulation/issues) page for this project on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
