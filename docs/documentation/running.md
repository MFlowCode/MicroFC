# Running

Python input file `input.py` defines dependencies and logistics, and input parameters for each simulation case.
In this section, details of the input file and how to edit it are described.
The user can also leverage the example input files as necessary.

## Parameters

There are multiple sets of parameters that must be specified in the python input file:
1. [Runtime Parameters](#1-runtime)
2. [Computational Domain Parameters](#2-computational-domain)
3. [Patch Parameters](#3-patches)
4. [Fluid Material's Parameters](#4-fluid-materials)
5. [Simulation Algorithm Parameters](#5-simulation-algorithm)
6. [Formatted Database and Structure Parameters](#6-formatted-output)

Items 7 and 8 are optional sets of parameters that activate the acoustic source model and ensemble-averaged bubble model, respectively.
Definition of the parameters is described in the following subsections.

### 1. Runtime

| Parameter        | Type           | Description                      |
| ---:             |    :----:      |          :---                    |
| `case_dir`       | String         | Case script directory            |
| `run_time_info`  | Logical        | Output run-time information      |

- `case_dir` specifies the directory where the python input file is located.

- `run_time_info` generates a text file that includes run-time information including the CFL number(s) at each time-step.

### 2. Computational Domain

| Parameter                | Type    | Description                      |
| ---:                     | :----:  |          :---                    |
| `x[y]_domain%beg[end]` | Real    | Beginning [ending] of the $x$[y]-direction domain    |
| `stretch_x[y]`         | Logical | Stretching of the mesh in the $x$[y]-direction |
| `a_x[y]`               | Real    | Rate at which the grid is stretched in the $x$[y]-direction |
| `x[y]_a`               | Real    | Beginning of the stretching in the negative $x$[y]-direction |
| `x[y]_b`               | Real    | Beginning of the stretching in the positive $x$[y]-direction |
| `m`                      | Integer | Number of grid cells in the $x$-coordinate direction |
| `n`                      | Integer | Number of grid cells in the $y$-coordinate direction |
| `dt`                     | Real    | Time step size |
| `t_step_start`           | Integer | Simulation starting time step |
| `t_step_stop`            | Integer | Simulation stopping time step |
| `t_step_save`            | Integer | Frequency to output data |

The parameters define the boundaries of the spatial and temporal domains, and their discritization that are used in simulation.

- `[x,y]_domain%[beg,end]` define the spatial domain in $x$, $y$, and $z$ Cartesian coordinates:

$$ x \in \left[ x \\_ domain \\% beg, x \\_ domain \\% end \right], y \in \left[ y \\_ domain \\% beg, y \\_ domain \\% end \right], z \in \left[ z \\_ domain \\% beg, z \\_ domain \\% end \right] $$

- $m$, $n$, and $p$ define the number of finite volume cells that uniformly discritize the domain along the $x$, $y$, and $z$ axes, respectively.
Note that the actual number of cells in each coordinate axis is given as $[m,n,p]+1$.
For example, $m=n=p=499$ discretizes the domain into $500^3$ cells. 
When the simulation is 2D/axi-symmetric or 1D, it requires that $p=0$ or $p=n=0$, respectively.

- `stretch_[x,y]` activates grid stretching in the $[x,y]$ directions.
The grid is gradually stretched such that the domain boundaries are pushed away from the origin along a specified axis.

- `a_[x,y]`, `[x,y]_a`, and `[x,y]_b` are parameters that define the grid stretching function. When grid stretching along the $x$ axis is considered, the stretching function is given as:

$$ x_{cb,stretch} = x_{cb} + \frac{x_{cb}}{a_x} \Bigg[ \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_{cb}-x_a)}{L} \right) \right] + \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_{cb}-x_b)}{L} \right) \right] -2 \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_b-x_a)}{2L} \right) \right]  \Bigg] $$

where `x_cb` and `x_[cb,stretch]` are the coordinates of a cell boundary at the original and stretched domains, respectively. `L` is the domain length along the `x` axis: `L`=`x_domain%end`-`x_domain%beg`. Crudely speaking, `x_a` and `x_b` define the coordinates at which the grid begins to get stretched in the negative and positive directions along the $x$ axis, respectively. $a_x$ defines the smoothness of the stretching. Stretching along the $y$ and $z$ axes follows the same logistics. Optimal choice of the parameters for grid stretching is case-dependent and left to the user.

- `dt` specifies the constant time step size that is used in simulation. The value of `dt` needs to be sufficiently small such that the Courant-Friedrichs-Lewy (CFL) condition is satisfied.

- `t_step_start` and `t_step_end` define the time steps at which simulation starts and ends, respectively. `t_step_save` is the time step interval for data output during simulation. To newly start simulation, set `t_step_start`=0. To restart simulation from $k$-th time step, set `t_step_start`=k.

### 3. Patches

| Parameter           | Type    | Description                                                  |
| ---:                | :----:  |          :---                                                |
| `num_patches`       | Integer | Number of initial condition geometric patches.               |
| `num_fluids`	      | Integer | Number of fluids/components present in the flow.             |
| `geometry` *        | Integer | Geometry configuration of the patch.                         |
| `alter_patch(i)` *  | Logical | Alter the $i$-th patch.                                      |
| `x[y]_centroid` * | Real    | Centroid of the applied geometry in the $[x,y]$-direction. |
| `length_x[y]` *   | Real    | Length, if applicable, in the $[x,y]$-direction.           |
| `radius` *          | Real	   | Radius, if applicable, of the applied geometry.              |
| `smoothen` *        | Logical | Smoothen the applied patch.                                  |
| `smooth_patch_id` * | Integer | A patch with which the applied patch is smoothened.          |
| `smooth_coeff` *    | Real    | Smoothen coefficient.                                        |
| `alpha(i)` *        | Real    | Volume fraction of fluid $i$.                                |
| `alpha_rho(i)` *    | Real    | Partial density of fluid $i$.                                |
| `pres` *            | Real    | Pressure.                                                    |
| `vel(i)` *          | Real    | Velocity in direction $i$.                                   |

*: These parameters should be prepended with `patch_icpp(j)%` where $j$ is the patch index. 

The Table lists the patch parameters. The parameters define the geometries and physical parameters of fluid components (patch) in the domain at initial condition. Note that the domain must be fully filled with patche(s). The code outputs error messages when an empty region is left in the domain.

- `num_patches` defines the total number of patches defined in the domain. The number has to be a positive integer.

- `num_fluids` defines the total number of fluids defined in each of the patches. The number has to be a positive integer.

- `patch_icpp(j)%geometry` defines the type of geometry of $j$-th patch by using an integer from 1 to 13. Definition of the patch type for each integer is listed in table [Patch Types](#patch-types).

- `[x,y]_centroid`, `length_[x,y]`, and/or `radius` are used to uniquely define the geometry of the patch with given type. Requisite combinations of the parameters for each type can be found in is listed in table [Patch types](#patch-types).

- `smoothen` activates smoothening of the boundary of the patch that alters the existing patch.
When smoothening occurs, fluids of the two patches are mixed in the region of the boundary.
For instance, in the aforementioned case of the cylindrical patch immersed in the rectangular patch, smoothening occurs when `patch_icpp(2)smoothen`=TRUE. `smooth_coeff` controls the thickness of the region of smoothening (sharpness of the mixture region). The default value of `smooth_coeff` is unity. The region of smoothening is thickened with decreasing the value.
Optimal choice of the value of `smooth_coeff` is case-dependent and left to the user.

- `patch_icpp(j)alpha(i)`, `patch_icpp(j)alpha_rho(i)`, `patch_icpp(j)pres`, and `texttt{patch_icpp(j)vel(i)` define for $j$-th patch the void fraction of `fluid(i)`, partial density of `fluid(i)`, the pressure, and the velocity in the $i$-th coordinate direction. These physical parameters must be consistent with fluid material's parameters defined in the next subsection.
See also `adv_alphan` in table [Simulation Algorithm Parameters](#5-simulation-algorithm).

### 4. Fluid Material’s

| Parameter | Type   | Description                                    |
| ---:      | :----: |          :---                                  |
| `gamma`   | Real   | Stiffened-gas parameter $\Gamma$ of fluid.     |
| `pi_inf`  | Real   | Stiffened-gas parameter $\Pi_\infty$ of fluid. |
| `Re(1)`   | Real   | Shear viscosity of fluid.                      |
| `Re(2)`   | Real   | Volume viscosity of fluid.                     |

Fluid material's parameters. All parameters should be prepended with `fluid_pp(i)` where $i$ is the fluid index.

The table lists the fluid material's parameters.
The parameters define material's property of compressible fluids that are used in simulation.

- `fluid_pp(i)%gamma` and `fluid_pp(i)%pi_inf` define $\Gamma$ and $\Pi$ as parameters of $i$-th fluid that are used in stiffened gas equation of state.

- `fluid_pp(i)%Re(1)` and `fluid_pp(i)%Re(2)` define the shear and volume viscosities of $i$-th fluid, respectively.
When these parameters are undefined, fluids are treated as inviscid.
Details of implementation of viscosity in MFC can be found in [Coralic (2015)](references.md#Coralic15).

### 5. Simulation Algorithm

| Parameter              | Type    | Description                                    |
| ---:                   | :----:  |          :---                                  |
| `bc_[x,y]\%beg[end]` | Integer | Beginning [ending] boundary condition in the $[x,y]$-direction (negative integer, see table [Boundary Conditions](#boundary-conditions)) |
| `time_stepper`         | Integer | Runge--Kutta order [1--5] |
| `weno_vars`	           | Integer | WENO reconstruction on [1] Conservative; [2] Primitive variables |
| `weno_order`	         | Integer | WENO order [1,3,5] |
| `weno_eps`	           | Real    | WENO perturbation (avoid division by zero) |

The table lists simulation algorithm parameters.
The parameters are used to specify options in algorithms that are used to integrate the governing equations of the multi-component flow based on the initial condition.
Models and assumptions that are used to formulate and discritize the governing equations are described in [Bryngelson et al. (2019)](references.md#Bryngelson19).
Details of the simulation algorithms and implementation of the WENO scheme can be found in [Coralic (2015)](references.md#Coralic15).

- `bc_[x,y]%[beg,end]` specifies the boundary conditions at the beginning and the end of domain boundaries in each coordinate direction by a negative integer from -1 through -12. See table [Boundary Conditions](#boundary-conditions) for details.

- `time_stepper` specifies the order of the Runge-Kutta (RK) time integration scheme that is used for temporal integration in simulation, from the 1st to 5th order by corresponding integer. 
Note that `time_stepper` $=$ 3 specifies the total variation diminishing (TVD), third order RK scheme ([Gottlieb and Shu, 1998](references.md#Gottlieb98)).

- `weno_vars` $=$ 1 and 2 correspond to conservative variables and primitive variables, respectively.

- `weno_order` specifies the order of WENO scheme that is used for spatial reconstruction of variables by an integer of 1, 3, and 5, that correspond to the 1st, 3rd, and 5th order, respectively.

- `weno_eps` specifies the lower bound of the WENO nonlinear weights. Practically, `weno_eps` $<10^{-6}$ is used.

### 6. Formatted Output

| Parameter            | Type    | Description                                    |
| ---:                 | :----:  |          :---                                  |
| `format`             | Integer | Output format. [1]: Silo-HDF5; [2] Binary	|
| `precision`          | Integer | [1] Single; [2] Double	 |
| `parallel_io`        | Logical | Parallel I/O	|
| `cons_vars_wrt`      | Logical | Write conservative variables \|
| `prim_vars_wrt`      | Logical | Write primitive variables	|
| `omega_wrt(i)`       | Logical | Add the $i$-direction vorticity to the database	 |
| `schlieren_wrt`      | Logical | Add the numerical schlieren to the database|
| `fd_order`           | Integer | Order of finite differences for computing the vorticity and the numerical Schlieren function [1,2,4] |
| `schlieren_alpha(i)` | Real    | Intensity of the numerical Schlieren computed via `alpha(i)` |
| `probe_wrt`          | Logical | Write the flow chosen probes data files for each time step	|
| `num_probes`         | Integer | Number of probes	|
| `probe(i)%[x,y]`   | Real	   | Coordinates of probe $i$	|

The table lists formatted database output parameters. The parameters define variables that are outputted from simulation and file types and formats of data as well as options for post-processing.

- `format` specifies the choice of the file format of data file outputted by MFC by an integer of 1 and 2. `format` $=$ 1 and 2 correspond to Silo-HDF5 format and binary format, respectively.

- `precision` specifies the choice of the floating-point format of the data file outputted by MFC by an integer of 1 and 2. `precision` $=$ 1 and 2 correspond to single-precision and double-precision formats, respectively.

- `parallel_io` activates parallel input/output (I/O) of data files. It is highly recommended to activate this option in a parallel environment.
With parallel I/O, MFC inputs and outputs a single file throughout pre-process, simulation, and post-process, regardless of the number of processors used.
Parallel I/O enables the use of different number of processors in each of the processes (i.e. simulation data generated using 1000 processors can be post-processed using a single processor).

- `cons_vars_wrt` and `prim_vars_wrt` activate output of conservative and primitive state variables into the database, respectively.

- `[variable's name]_wrt` activates output of the each specified variable into the database.

- `schlieren_alpha(i)` specifies the intensity of the numerical Schlieren of $i$-th component.

- `fd_order` specifies the order of finite difference scheme that is used to compute the vorticity from the velocity field and the numerical schlieren from the density field by an integer of 1, 2, and 4. `fd_order` $=$ 1, 2, and 4 correspond to the first, second, and fourth order finite difference schemes, respectively.

- `probe_wrt` activates output of state variables at coordinates specified by `probe(i)%[x;y]`.


## Enumerations

### Boundary conditions

| #    | Description | 
| ---: | :---        |
|  -1  | Periodic |
|  -2  | Reflective |
|  -3  | Ghost cell extrapolation |
|  -4  | Riemann extrapolation |
|  -5  | Slip wall |
	
The boundary condition supported by the MFC are listed in table [Boundary Conditions](#boundary-conditions). Their number (`#`)
corresponds to the input value in `input.py` labeled `bc_[x,y]%[beg,end]` (see table [Simulation Algorithm Parameters](#5-simulation-algorithm)).
The entries labeled "Characteristic." are characteristic boundary conditions based on [Thompson (1987)](references.md#Thompson87) and [Thompson (1990)](references.md#Thompson90).

### Patch types

| #    | Name           | Dim.  | Smooth | Description |
| ---: | :----:         | :---: |  :---: | :--- |
| 1    | Line segment 	| 1     | N      | Requires `x_centroid` and `x_length`. |
| 2    | Circle 		    | 2     | Y      | Requires `[x,y]_centroid` and `radius`. |
| 3    | Rectangle 	    | 2     | N      | Coordinate-aligned. Requires `[x,y]_centroid` and `[x,y]_length`. |
| 4    | Sweep line 		| 2     | Y      | Not coordinate aligned. Requires `[x,y]_centroid` and `normal(i)`. |
| 5    | Ellipse 		    | 2     | Y      | Requires `[x,y]_centroid` and `radii(i)`. |
| 6    | Vortex 		    | 2     | N      | Isentropic flow disturbance. Requires `[x,y]_centroid` and `radius`. |
| 7    | 2D analytical 	| 2     | N      | Assigns the primitive variables as analytical functions. |

The patch types supported by the MFC are listed in table [Patch Types](#patch-types). This includes
types exclusive to one-, two-, and three-dimensional problems. The patch type number (`#`)
corresponds to the input value in `input.py` labeled  `patch_icpp(j)%geometry` where
$j$ is the patch index. Each patch requires a different set of parameters, which are 
also listed in this table.

## Running

MFC can be run using `mfc.sh`'s `run` command. It supports both interactive and
batch execution, the latter being designed for multi-socket systems, namely supercomputers,
equipped with a scheduler such as PBS, SLURM, and LSF. A full (and updated) list
of available arguments can be acquired with `./mfc.sh run -h`. Example Python input
files can be found in the [examples/](examples/) directory. They print a Python dictionary containing input parameters for MFC. Their contents, and a guide to filling them out, are documented
in the user manual. A commented, tutorial script
can also be found in [examples/3d_sphbubcollapse/](examples/3D_sphbubcollapse/case.py).

The skeleton for an input file may look like the following:

```python
#!/usr/bin/env python3

import json

# Configuring case dictionary
print(json.dumps({
  # Insert case parameters here
  ...
}))
```

Thus, you can run your case file with Python to view the computed case dictionary
that will be processed by MFC when you run. This is particularly useful when
computations are done in Python to generate the case.

### Interactive Execution

To run all stages of MFC, that is [pre_process](src/pre_process/), [simulation](src/simulation/), and [post_process](src/post_process/)   on the sample case [2D_shockbubble](examples/2D_shockbubble/),

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py
```

If you want to run a subset of the available stages, you can use the `-t` argument.
To use multiple threads, use the `-n` option along with the number of threads you wish to use.
If a (re)build is required, it will be done automatically, with the number of threads
specified with the `-j` option.

For example,

- Running [pre_process](src/pre_process/) with 2 cores:

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -t pre_process -n 2
```

- Running [simulation](src/simulation/) and [post_process](src/post_process/)
using 4 cores:

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -t simulation post_process -n 4
```

Most parameters have sensible defaults which can be overridden in [defaults.yaml](defaults.yaml):

https://github.com/MFlowCode/MFC-develop/blob/d74e714b08562a9f8f815112e05df54c99c8c18f/defaults.yaml#L12-L21

On some computer clusters, MFC might select the wrong MPI program to execute your application
because it uses a general heuristic for its selection. Notably, `srun` is known to fail on some SLURM
systems when using GPUs or MPI implementations from different vendors, whereas `mpirun` functions properly. To override and manually specify which
MPI program you wish to run your application with, please use the `-b <program name>` option (i.e `--binary`).

Additional flags can be appended to the MPI executable call using the `-f` (i.e `--flags`) option.

Please refer to `./mfc.sh run -h` for a complete list of arguments and options, along with their defaults.

### Batch Execution

The MFC detects which scheduler your system is using and handles the creation and
execution of batch scripts. The batch engine is requested with the `-e batch` option.
Whereas the interactive engine can execute all of MFC's codes in succession, the batch engine
requires you to only specify one target with the `-t` option. The number of nodes and GPUs can, 
respectively be specified with the `-N` (i.e `--nodes`) and `-g` (i.e `--gpus-per-node`) options.

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -e batch -N 2 -n 4 -g 4 -t simulation
```

Other useful arguments include:

- `-# <job name>` to name your job. (i.e `--name`)
- `-@ sample@example.com` to receive emails from the scheduler. (i.e `--email`)
- `-w hh:mm:ss` to specify the job's maximum allowed walltime. (i.e `--walltime`)
- `-a <account name>` to identify the account to be charged for the job. (i.e `--account`)
- `-p <partition name>` to select the job's partition. (i.e `--partition`)

Since some schedulers don't have a standardized syntax to request certain resources, MFC can only
provide support for a restricted subset of common configuration options. If MFC fails
to execute on your system, or if you wish to adjust how the program runs and resources
are requested to be allocated, you are invited to modify the template batch script for your queue system.
Upon execution of `./mfc.sh run`, MFC fills in the template with runtime parameters, to
generate the batch file it will submit. These files are located in the [templates](templates/)
directory. To request GPUs, modification of the template will be required on most systems.

- Lines that begin with `#>` are ignored and won't figure in the final batch
script, not even as a comment.

- Statements of the form `${expression}` are string-replaced to provide
runtime parameters, most notably execution options. They reference the variables in the
same format as those under the "run" section of [defaults.yaml](defaults.yaml), replacing
`-` for `_`. You can perform therein any Python operation recognized by the built-in `expr()` function.

As an example, one might request GPUs on a SLURM system using the following:

```
#SBATCH --gpus=v100-32:{gpus_per_node*nodes}
```

- Statements of the form `{MFC::expression}` tell MFC where to place the common code,
across all batch files, that is required for proper execution. They are not intended to be
modified by users.

**Disclaimer**: IBM's JSRUN on LSF-managed computers does not use the traditional node-based approach to
allocate resources. Therefore, the MFC constructs equivalent resource-sets in task and GPU count.

### Example Runs

- Oak Ridge National Laboratory's [Summit](https://www.olcf.ornl.gov/summit/):

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -e batch    \
               -N 2 -n 4 -g 4 -t simulation -a <redacted>
```

- University of California, San Diego's [Expanse](https://www.sdsc.edu/services/hpc/expanse/):

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -e batch -p GPU -t simulation \
               -N 2 -n 8 -g 8 -f="--gpus=v100-32:16" -b mpirun –w 00:30:00
```
