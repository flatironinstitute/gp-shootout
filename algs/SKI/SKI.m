function [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrg, opts)
% SKI   GP regression via SKI algorithm in 1 or more dimension
%
% [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrg, opts)
%  performs Gaussian process regression via SKI using gpytorch.
%
%  NOTE: when performing regression for data on [0, 1], we need data points
%  (x) at x=0 and x=1.   *** is this still true?
%  
% Inputs:
%  x    - points (ordinates) where observations taken, d*N real array for d dims
%  meas - observations at the data points, length-N real (or complex?) array
%  sigmasq - noise variance at data points, nonnegative scalar
%  ker  - covariance kernel info with at least the field:
%         khat - kernel Fourier transform function handle, must act elementwise
%                on array of wavevector magnitudes |xi| (isotropic).
%  xtrg - [optional, or may be empty] targets points, d*n real array for d dims.
%         If non-empty, attempts to compute ytrg outputs.
%  opts - [optional] struct controlling method params including:
%         grid_size - ski grid size
%         default_grid_size - use ski default grid size
%
% Outputs:
%  y - struct with fields of regression results corresp to given data points x:
%     mean - posterior mean vector, N*1
%  ytrg - [optional; otherwise empty] struct of regression at new targets xtrg:
%     mean - posterior mean vector, n*1
%  info - diagnostic struct containing fields:
%     cputime - struct with timing fields including: total
%
% If called without arguments, does a self-test.

% Notes: 1) this code is a wrapper to a python implementation of SKI
if nargin==0, test_SKI; return; end
if nargin<5, xtrg = []; end
do_trg = ~isempty(xtrg);
if nargin<6, opts = []; end

[dim, N] = size(x);
if numel(meas)~=N, error('sizes of meas and x must match!'); end
% use default grid size if specified by user or if opts.grid_size doesn't
% exist
do_default = false;
if isfield(opts, 'default_grid_size')
    if opts.default_grid_size
        do_default = true;
    end
end
if ~isfield(opts,'grid_size'), do_default = true; opts.grid_size = 0; end     % default


xsol = [x, xtrg];  % hack for now which adds meas pts to target list
xpy = py.numpy.array(x');
ypy = py.numpy.array(meas');
testxpy = py.numpy.array(xsol');
do_double = true;
ski_out = py.ski.gpr(xpy, ypy, testxpy, opts.grid_size, sigmasq, ker.fam, ker.l, do_double, do_default);
% unpack ski output
yhat = double(ski_out{1})';
info.cpu_time.total = ski_out{2};

y.mean = yhat(1:N);   % hack for now to split out posterior means into two types
ytrg.mean = yhat(N+1:end);


%%%%%%%%%%
function test_SKI

% data
dim = 1;
N = 10;
x = rand(dim, N);
x(1) = 0;
x(2) = 1;
sigma_true = 1.0;
meas = cos(x) + sigma_true * rand(1, N);

% test points
ntest = 10;
testx = sort(rand(dim, ntest));
testx = linspace(0, 1, ntest);
testx(1) = x(1);
testx(ntest) = x(N);

opts.grid_size = 100;
l = 0.2;
ker = SE_ker(dim, l);
sigmasq = sigma_true^2;
[yhat, ytrg, info] = SKI(x, meas, sigmasq, ker, testx, opts);
[yhat2, ytrg2, ~] = naive_gp(x, meas, sigmasq, ker, testx, []);
fprintf('SKI test dim=%d, N=%d, ntest=%d:   \tmax diff %.5g   \t(time %.3g s)\n', dim, N, ntest, max(abs(ytrg.mean - ytrg2.mean)), info.cpu_time.total);

ntest = 10000;
testx = linspace(0, 1, ntest);
testx(1) = x(1);
testx(ntest) = x(N);
[yhat, ytrg, info] = SKI(x, meas, sigmasq, ker, testx, opts);
[yhat2, ytrg2, ~] = naive_gp(x, meas, sigmasq, ker, testx, []);
fprintf('SKI test dim=%d, N=%d, ntest=%d:   \tmax diff %.5g   \t(time %.3g s)\n', dim, N, ntest, max(abs(ytrg.mean - ytrg2.mean)), info.cpu_time.total);


% now in 2d
% data
dim = 2;
N = 10;
x = rand(dim, N);
x(:, 1) = 0;
x(:, 2) = 1;
sigma_true = 1.0;
meas = cos(x(1, :) + 3 * x(2, :)) + sigma_true * rand(1, N);

% test points
ntest = 10;
testx = sort(rand(dim, ntest));

opts.grid_size = 100;
l = 0.2;
ker = SE_ker(dim, l);
sigmasq = sigma_true^2;
[yhat, ytrg, info] = SKI(x, meas, sigmasq, ker, testx, opts);
[yhat2, ytrg2, ~] = naive_gp(x, meas, sigmasq, ker, testx, []);
fprintf('SKI test dim=%d, N=%d, ntest=%d:   \tmax diff %.5g   \t(time %.3g s)\n', dim, N, ntest, max(abs(ytrg.mean - ytrg2.mean)), info.cpu_time.total);


% generates warnings: /usr/local/opt/python/Frameworks/Python.framework/
% Versions/3.7/lib/python3.7/site-packages/gpytorch/utils/interpolation.py:
% 119: UserWarning: Invalid MKL_NUM_THREADS variable value, stoi: no 
% conversion (Triggered internally at  ../aten/src/ATen/ParallelCommon.cpp
% :38.) left_boundary_pts = (lower_grid_pt_idxs < 0).nonzero(as_tuple=False)
