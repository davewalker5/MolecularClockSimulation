[![GitHub issues](https://img.shields.io/github/issues/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/issues)
[![Releases](https://img.shields.io/github/v/release/davewalker5/MolecularClockSimulation.svg?include_prereleases)](https://github.com/davewalker5/MolecularClockSimulation/releases)
[![License](https://img.shields.io/badge/License-mit-blue.svg)](https://github.com/davewalker5/MolecularClockSimulation/blob/main/LICENSE)
[![Language](https://img.shields.io/badge/language-python-blue.svg)](https://www.python.org)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/davewalker5/MolecularClockSimulation)](https://github.com/davewalker5/MolecularClockSimulation/)

# Molecular Clock Simulation

<img src="https://github.com/davewalker5/MolecularClockSimulation/blob/main/diagrams/relaxed-clock-explorer-simulation.png" width="100%">

A Python toolkit for simulating sequence evolution, analysing evolutionary distances and exploring distance-based phylogenetic reconstruction from first principles.

## Overview

Molecular clocks use DNA sequence variation to estimate evolutionary relationships and divergence times between species.

This project explores the computational foundations of molecular clock analysis by implementing the core algorithms from first principles. It combines reusable simulation engines, analytical tools and interactive visualisations to provide an educational framework for understanding sequence evolution and phylogenetic reconstruction.

Detailed documentation, implementation notes and usage guides are available in the project Wiki.

## Current Features

| Component                       | Purpose                                                            |
| ------------------------------- | ------------------------------------------------------------------ |
| Strict Molecular Clock Engine   | Simulate ultrametric evolutionary histories                        |
| Relaxed Molecular Clock Engine  | Simulate lineage-specific rate variation                           |
| Distance Matrix Calculator      | Estimate evolutionary distances using multiple substitution models |
| UPGMA Reconstruction            | Reconstruct rooted molecular-clock trees from distance matrices    |
| Neighbor Joining Reconstruction | Reconstruct distance-based trees without a strict-clock assumption |
| Interactive Explorers           | Combine the complete workflow into an educational environment      |

## Documentation

Further documentation, including the documentation on intended field usage and how to run the application, is available in the project [Wiki](https://github.com/davewalker5/MolecularClockSimulation/wiki).

## Feedback

To file issues or suggestions, please use the [Issues](https://github.com/davewalker5/MolecularClockSimulation/issues) page for this project on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
