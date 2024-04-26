[![status](https://github.com/olegrog/kesolver/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/olegrog/kesolver/actions/workflows/ubuntu.yml)

Kinetic Equation Solver (kesolver) is a solver of the Boltzmann equation
using the conservative projection method and its further extensions.

Tcheremissine, F. G. "Conservative evaluation of Boltzmann collision integral in discrete
ordinates approximation." Computers & Mathematics with Applications 35.1-2 (1998): 215-221.

## Dependencies

Mandatory:

1. A C++ compiler supporting C++11 standard.
2. [JsonCpp](https://github.com/open-source-parsers/jsoncpp) library.
See instructions [here](JsonCpp/README.md).

Optionally:

1. An MPI wrapper for the C++ compiler.
2. [Gmsh](https://gmsh.info/) for mesh generation.
3. [ParaView](https://www.paraview.org/) for data visualization.

## Building

Run this command to create 2 executables: `kes` and `kes-mpi`.

```shell
cd src && make -j
```

By default, `mpic++` is used to build an MPI version of the solver.

Currently, only the x86-64 architecture is supported.
For a Mac with an arm-based processor, use `arch -x86_64 make -j`.

## Running

The executable files expect a config file in a special `kei` format,
which includes a prepared computational mesh.
The path to it is given as a single argument, e.g.

```shell
./kes config.kei
```

A parallel version can be run as follows:

```shell
mpirun ./kes-mpi config.kei
```

Examples of configuration files can be found in [tutorials](tutorial).
