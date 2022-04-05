function [yhat, ytrg, info] = naive_gp(x, meas, sigmasq, ker, xtrg, opts)
% NAIVE_GP   Perform naive slow O(N^3) GP regression, arbitrary dimension.
%
% [y, ytrg, info] = naive_gp(x, meas, sigmasq, ker, xtrg, opts)
%  performs Gaussian process regression using an arbitrary
%  translation-invariant isotropic prior kernel ker, conditioned on
%  the y-values meas at the set of points x, in arbitrary spatial
%  dimension. A naive O(N^3) algorithm is used, without any
%  approximations apart from rounding error. It is thus a reference
%  implemention.
%
% Inputs:
%  x    - points (ordinates) where observations taken, d*N real array for d dims
%  meas - observations at the data points, length-N real (or complex?) array
%  sigmasq - noise variance at data points, nonnegative scalar
%  ker  - covariance kernel info with at least the field:
%         k - kernel function handle, must act elementwise on array of distances
%  xtrg - [optional, or may be empty] targets points, d*n real array for d dims.
%         If non-empty, attempts to compute ytrg outputs.
%  opts - [optional] struct controlling method params including:
%         getcovar - if true, compute conditional covariances
%
% Outputs:
%  y - struct with fields of regression results corresp to given data points x:
%     mean - posterior mean vector, N*1
%     alpha - N*1 vector of weights for reconstruction of mean
%     covar - posterior covariance, N*N. [optional]
%  ytrg - [optional; otherwise empty] struct of regression at new targets xtrg:
%     mean - posterior mean vector, n*1, computed stably
%     meanbook - posterior mean vector, n*1, computed by K*alpha
%     covar - posterior covariance, n*n. [optional]
%  info - diagnostic struct containing fields:
%     cputime - list of times in seconds for fill, solve + y eval, trg eval, etc.
%
% If called without arguments, does a self-test
if nargin==0, test_naive_gp; return; end
do_trg = (nargin>=5 && ~isempty(xtrg));
[dim,N] = size(x);
meas = meas(:);
if numel(meas)~=N, error('sizes of meas and x must match!'); end
if N>1e4, warning('N getting too big for naive method!'); end

tic; K = densekermat(ker.k,x);
info.cputime(1) = toc;

tic;
alpha = (K + sigmasq*eye(N)) \ meas;       % dense solve, probably
y.mean = meas - sigmasq*alpha;             % magic stable cheap formula
y.meanbook = K*alpha;                      % the book formula, slow, unstable
info.cputime(2) = toc;

if do_trg
  [~,n] = size(xtrg);        % # targs
  tic;
  B = densekermat(ker.k,xtrg,x);   %  n-by-N
  yhat = B * alpha;    % do the GP sum, naively (not as stable as above)
  info.cputime(3) = toc;
end


%%%%%%%%%%
function test_naive_gp   % basic tests for now
N = 1e3;        % problem size
l = 0.1;        % SE kernel scale
ker.k = @(d) exp(-(0.5/l^2)*(d.*d));
sigma = 0.01;   % for regression
sigmadata = sigma;
freqdata = 3;   % how oscillatory underlying func; freq >> 0.3/l misspecified

for dim = 1:2
  fprintf('test naive_gp, dim=%d...\n',dim)
  x = rand(dim,N);     % in [0,1]^dim
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;
  f = @(x) cos(2*pi*wavevec*x + 1.3);   % underlying func, will give row vec
  meas = f(x) + sigmadata*randn(size(x));   % noisy data
  [y, ~, info] = naive_gp(x, meas, sigma^2, ker.k);
  fprintf('CPU times:'); disp(info.cputime)
  fprintf('y.mean: rms resid of lin sys = %.3g\n', norm(y.mean-y.meanbook)/sqrt(N))
end
