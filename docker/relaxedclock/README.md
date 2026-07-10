# relaxedclockexplorer

The Relaxed Molecular Clock Explorer provides an interactive environment for exploring the relaxed molecular clock simulation from the [_Molecular Clock Simulation_](https://github.com/davewalker5/MolecularClockSimulation) project and the complete workflow from sequence evolution to phylogenetic reconstruction.

Built on top of the project’s Relaxed Molecular Clock Engine, the explorer allows users to generate evolutionary histories, analyse simulated DNA sequences using multiple evolutionary distance models, visualise distance matrices and reconstruct phylogenetic trees using UPGMA or Neighbor Joining.

This image allows the explorer to be run without cloning or building the project locally, making it a convenient way to explore the simulation and reconstruction workflow.

## Getting Started

### Prerequisites

In order to run this image you'll need docker installed.

- [Windows](https://docs.docker.com/windows/started)
- [OS X](https://docs.docker.com/mac/started/)
- [Linux](https://docs.docker.com/linux/started/)

### Usage

#### Container Parameters

The following "docker run" parameters are recommended when running the relaxedclockexplorer image:

| Parameter  | Value        | Purpose                                                 |
| ---------- | ------------ | ------------------------------------------------------- |
| -d         | -            | Run as a background process                             |
| -p         | 80:5000      | Expose the container's port 5000 as port 80 on the host |
| --platform | linux/amd64  | Target architecture ; this must be linux/amd64          |
| --name     | relaxedclock | Name of the container once running                      |
| --rm       | -            | Remove the container automatically when it stops        |

For example:

```shell
docker run -d -p 80:8501 --platform linux/amd64 --name relaxedclock --rm davewalker5/relaxedclockexplorer:latest
```

The port number "80" can be replaced with any available port on the host.

#### Running the Image and Accessing the Application

To run the image, enter the following commands, substituting an available port on the local machine for port 80, as described:

```shell
docker run -d -p 80:8501 --platform linux/amd64 --name relaxedclock --rm davewalker5/relaxedclockexplorer:latest
```

Once the container is running, browse to the following URL on the host:

http://localhost:80

## Project

Source code, documentation and release notes are available on GitHub:

- [MolecularClockSimulation on GitHub](https://github.com/davewalker5/MolecularClockSimulation)

## Versioning

For the versions available, see the [tags on this repository](https://github.com/davewalker5/MolecularClockSimulation/tags).

## Authors

- **Dave Walker** - _Initial work_ -

See also the list of [contributors](https://github.com/davewalker5/MolecularClockSimulation/contributors) who
participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/davewalker5/MolecularClockSimulation/blob/master/LICENSE) file for details.
