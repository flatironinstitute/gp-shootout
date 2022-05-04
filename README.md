# gp-shootout

Benchmark and compare several large-scale GP regression methods in 1D, 2D, and 3D,
including our implementation of the equispaced Fourier method (EFGP).
We also generate figures and tables for the paper.

Authors: Philip R Greengard, Alex H Barnett, Manas Rachh.


### Installation

`git clone --recurse-submodules https://github.com/flatironinstitute/gp-shootout.git`

This will also download install some submodule packages (currently: RLCM, FLAM).
In addition your system must also have the following.
Required dependencies:

* MATLAB (tested on R2021b)

Dependencies specific to methods:

* For EFGP (and its option SSGP): FINUFFT (version 2.0 or later; please specify its location in `startup.m`).
* For SKI: Python 3.7, 3.8, or 3.9 (not 3.10, since MATLAB has not caught up), with
python environment in which MATLAB was opened to include:
   - numpy
   - torch
   - gpytorch
   - pytorch
* For RLCM: C++ compiler (and RLCM installed as submodule)
* FLAM: (FLAM is installed as a submodule)

To test the basic installation, start MATLAB from the top-level `gp-shootout`
directory (which will execute `startup`), then within MATLAB type `test_all`.

Advanced: to build then test all wrapped non-MATLAB methods:

1) make sure you can call python from matlab, eg via `py.sys.version`
2) from shell do `(cd algs/RLCM; ./buildit.sh)`

then from MATLAB run `test_all_nonmatlab`.



### Usage

If you did not start MATLAB from the top-level directory, then run `startup` to add required paths and apply useful settings.

Look in `drivers` for example scripts. You may try to run `expt` for a demo.

to do: minimally complete example...


### To do

* add Matern to RLCM (it's SE only so far)
* FLAMGP in 3D
* insert switch from EFGP to SSGP w/ same xi quad nodes.
* explore EFGP accel via padding fftn to powers of only 2,3,5.
* other top-level fig-generating driver scripts
* understand empirical error breakdown as sigma->0 at large N
* understand CG num iters growing like 1/sqrt(tol)
* protect from annoying FLAM chol fail if tol too large


### Done (CHANGELOG)

* rationalize interface for Philip dim-specific codes, move `getL()` into `kernels`
* add Matern nu=3/2, 5/2 in all 3 dims and to tester - a Matern kernel with nu as an input.
* CO2 dataset, 2D, N=1e6
* GPytorch (SKI) implementation
* Heaton et al. 2019 comparison datasets (two of), 2D, N=1e5
* posterior variances at target points for gp_naive
* testing routines for comparing methods to gp_naive for small problems
* rationalized dims 1,2,3 EFGP, handle arb shifted & scaled x data
* covar output to `naive_gp`?
* datasets -> `data/*`
* other methods -> `algs/*` (Python via system calls from matlab)
* 1d option for dense linear solve (vs. iterative)
* RLCM wrapped (SE ker only) via binary tmp file IO, tested in d=1,2,3.
* switched matlab pyenv engine setting to OutOfProcess to prevent crashes
* SKI fixed mem alloc problem
