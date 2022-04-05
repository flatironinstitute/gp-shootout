# gp-shootout

Benchmark and compare large-scale GP regression methods in 1D, 2D, and 3D,
including our implementation of the equispaced Fourier method (EFGP).
We also generate figures and tables for the paper.

Authors: Phillip R Greengard, Alex H Barnett, Manas Rachh.

### Installation

Required dependencies:

* MATLAB (tested on R2021b)

Method dependencies:

* For EFGP: FINUFFT (version 2.0 or later)
* For GPyTorch: Python with PyTorch and GPyTorch
* ...

To test your installation from MATLAB shell: `startup; test_all`

### Usage

From MATLAB, run `startup` to add required paths and apply useful settings.


