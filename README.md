# gp-shootout

Benchmark and compare large-scale GP regression methods in 1D, 2D, and 3D,
including our implementation of the equispaced Fourier method (EFGP).
We also generate figures and tables for the paper.

Authors: Philip R Greengard, Alex H Barnett, Manas Rachh.


### Installation

`git clone --recurse-submodules https://github.com/flatironinstitute/gp-shootout.git`

Your system must also have the following.
Required dependencies:

* MATLAB (tested on R2021b)

Method dependencies:

* For EFGP: FINUFFT (version 2.0 or later; please specify its location in `startup.m`)
* For GPyTorch: Python with PyTorch and GPyTorch
* ...

To test your installation from MATLAB shell: `startup; test_all`


### Usage

From MATLAB, first run `startup` to add required paths and apply useful settings.

to do: Minimally complete example...


### To do

* make error benchmarks for methods, vs naive_gp, for N<=1e4 probs
   ... or just use dense K mat to check lin sys residual for methods.
* other top-level fig-generating driver scripts
* ~~rationalize interface for Philip dim-specific codes, move `getL()` into `kernels`~~
* add Matern nu=3/2, 5/2 in all 3 dims and to tester (PG: I think we'll want to have a Matern kernel with nu as an input) 
* add covar and new-targ output to `naive_gp`
* datasets -> `data/*`.   Eg Heaton competition from 2019, 2D, N~1e5.
* other methods -> `algs/*` (Python via system calls from matlab?)

