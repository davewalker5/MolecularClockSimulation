"""Compatibility exports for the Molecular Clock Simulation package."""

from strictclock import (
    SimulationConfig,
    SimulationResult,
    TreeNode,
    load_config,
    run_simulation,
    write_outputs,
)

__all__ = [
    "SimulationConfig",
    "SimulationResult",
    "TreeNode",
    "load_config",
    "run_simulation",
    "write_outputs",
]
