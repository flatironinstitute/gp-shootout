# gp-shootout

Benchmark and compare large-scale GP regression methods in 1D, 2D, and 3D,
including our implementation of the equispaced Fourier method (EFGP).
We also generate figures and tables for the paper.

Authors: Philip R Greengard, Alex H Barnett, Manas Rachh.


### Installation

Required dependencies:

* MATLAB (tested on R2021b)

Method dependencies:

* For EFGP: FINUFFT (version 2.0 or later)
* For GPyTorch: Python with PyTorch and GPyTorch
* ...

First adjust locations of MATLAB interfaces to external libs in `startup.m`.

Then, to test your installation from MATLAB shell: `startup; test_all`


### Usage

From MATLAB, first run `startup` to add required paths and apply useful settings.

to do: Minimally complete example...


### To do

* bring in Philip codes (wrap them for one interface) to algs/
* add Matern nu=3/2, 5/2 in all 3 dims and to tester
* add covar and new-targ output to naive_gp

