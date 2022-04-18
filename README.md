# gp-shootout

Benchmark and compare large-scale GP regression methods in 1D, 2D, and 3D,
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

* For EFGP: FINUFFT (version 2.0 or later; please specify its location in `startup.m`).
* For SKI: Python 3.7, 3.8, or 3.9 (not 3.10, since MATLAB has not caught up), with
python environment in which MATLAB was opened to include:
   - numpy
   - torch
   - gpytorch
   - pytorch
* For RLCF: C++ compiler, ...
* ...

To test your installation from MATLAB shell: `startup; test_all`


### Usage

From MATLAB, first run `startup` to add required paths and apply useful settings.

You may try to run `expt` for a demo.

to do: minimally complete example...


### To do

* make error benchmarks for methods, vs naive_gp, for N<=1e4 probs
   ... or just use dense K mat to check lin sys residual for methods.
* other top-level fig-generating driver scripts
* covar output to `naive_gp`?
* datasets -> `data/*`
* other methods -> `algs/*` (Python via system calls from matlab?)


### Done (CHANGELOG)

* rationalize interface for Philip dim-specific codes, move `getL()` into `kernels`
* add Matern nu=3/2, 5/2 in all 3 dims and to tester - a Matern kernel with nu as an input.
* CO2 dataset, 2D, N=1e6
* GPytorch (SKI) implementation
* Heaton et al. 2019 comparison datasets (two of), 2D, N=1e5
* posterior variances at target points for gp_naive
