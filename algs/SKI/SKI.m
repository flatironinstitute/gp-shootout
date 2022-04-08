function [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrg, opts)
% SKI   GP regression via SKI algorithm in dim=1,2
%
% [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrg, opts)
%  performs Gaussian process regression 
%  NOTE: when performing regression for data on [0, 1], we need data points
%  at x=0 and x=1
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
%
% Outputs:
%  y - struct with fields of regression results corresp to given data points x:
%     mean - posterior mean vector, N*1
%  ytrg - [optional; otherwise empty] struct of regression at new targets xtrg:
%     mean - posterior mean vector, n*1
%  info - diagnostic struct containing fields:
%     xis - Fourier xi nodes use
%     beta - m*1 vector of weight-space (Fourier basis) weights
%     cputime - list of times in seconds for various steps
%     iter - # iterations needed
%
% If called without arguments, does a self-test.


% Notes: 1) this code is a wrapper to a python implementation of ski
if nargin==0, test_SKI; return; end
if nargin<5, xtrg = []; end
do_trg = ~isempty(xtrg);
if nargin<6, opts = []; end

[dim, N] = size(x);
if ~isfield(opts,'grid_size'), opts.grid_size = N; end     % default
if numel(meas)~=N, error('sizes of meas and x must match!'); end

xsol = [x, xtrg];  % hack for now which adds meas pts to target list
                    % and transpose to Philip n*d shape
xpy = py.numpy.array(x);
ypy = py.numpy.array(meas);
testxpy = py.numpy.array(xsol);
ski_out = py.ski.gpr(xpy, ypy, testxpy, opts.grid_size, sigmasq, ker.fam, ker.l);
% unpack ski output
yhat = double(ski_out{1})';
info.cpu_time = ski_out{2};

y.mean = yhat(1:N);   % hack for now to split out posterior means into two types
ytrg.mean = yhat(N+1:end);


%%%%%%%%%%
function test_SKI

% data
dim = 1;
N = 10;
x = linspace(0, 1, N);
x = rand(N, dim)';
x(1) = 0;
x(N) = 1;
x = sort(x);
%N = 2;
%x = [0, 1];
meas = cos(x) + rand(1, N);

% test points
ntest = 10;
testx = sort(rand(dim, ntest));
testx(1) = x(1);
testx(ntest) = x(N);

opts.grid_size = 100;
l = 0.2;
ker = SE_ker(1, l);
sigmasq = 1.0;
[yhat, ytrg, info] = SKI(x, meas, sigmasq, ker, testx, opts);

% sigmasq = 1.0;
% % kern_family must be one of 'matern12', 'matern32', 'squared-exponential'
% kern_family = 'matern12';
% kern_family = 'squared-exponential';
% xpy = py.numpy.array(x);
% ypy = py.numpy.array(y);
% testxpy = py.numpy.array(testx);
% out1 = py.ski.gpr(xpy, ypy, testxpy, grid_size, sigma2, kern_family, l);

% compare to naive approach
ker = SE_ker(1, l);
[yhat2, ytrg2, info] = naive_gp(x, meas, sigmasq, ker, testx, []);

fprintf('1d max difference at target points %g \n', max(abs(ytrg.mean - ytrg2.mean)));
