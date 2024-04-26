# 1-D heat-transfer problem

This tutorial describes a solution of the Boltzmann equation
for a monatomic gas confined between parallel plates with constant temperatures
$T_1=0.9$ and $T_2=1.1$. The initial distribution function is the Maxwellian
with zero velocity, unit density and temperature.

Run this command

```shell
./run_all.sh
```

to generate a mesh using Gmsh, to run a solver, and convert the results to VTK format.
The latter can be visualized in ParaView.

| File | Description |
| --- | --- |
| `heat.geo` | Text file written in the Gmsh scripting language for mesh generation |
| `heat.kep_with_comments` | Commented configuration file written in JSON format |
| `run_all.sh` | Shell script used to run the computational pipeline |
| `clean_case.sh` | Shell script used to remove all the result and temporary files |
