function [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrg, opts)
% EFGP   GP regression via equispaced Fourier iterative method, in dim=1,2 or 3
%
% [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrg, opts)
%  performs Gaussian process regression using an arbitrary
%  translation-invariant isotropic prior kernel ker (via its Fourier transform
%  ker.khat), conditioned on the y-values meas at the set of points x, in
%  spatial dimension 1,2 or 3. The method is efficient and accurate only when
%  khat decays rapidly in Fourier space.
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
%         tol - desired tolerance, eg 1e-6
%
% Outputs:
%  y - struct with fields of regression results corresp to given data points x:
%     mean - posterior mean vector, N*1
%     beta - m*1 vector of weight-space (Fourier basis) weights
%  ytrg - [optional; otherwise empty] struct of regression at new targets xtrg:
%     mean - posterior mean vector, n*1
%  info - diagnostic struct containing fields:
%     xis - Fourier xi nodes used
%     cputime - list of times in seconds for various steps
%     iter - # iterations needed
%
% If called without arguments, does a self-test

% Note: this code is a wrapper to separate dimension functions
if nargin==0, test_EFGP; return; end
do_trg = (nargin>=5 && ~isempty(xtrg));
[dim,N] = size(x);
if numel(meas)~=N, error('sizes of meas and x must match!'); end

if dim==1
  [y.beta, info.xis, y.mean, info.iter, time_info] = function_space1d(x, y, sigma2, kern, eps, xsol)
elseif dim==2

elseif dim==3

else
  error('dim must be 1,2, or 3!');
end


%%%%%%%%%%
function test_EFGP   % basic tests for now, duplicates naive_gp *** to unify
N = 3e3;        % problem size
l = 0.1;        % SE kernel scale
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-10;

for dim = 1:2   % ..........
  fprintf('\ntest EFGP, sigma=%.3g, tol=%.3g, dim=%d...\n',sigma,opts.tol,dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  ker = SE_ker(dim,l);
  [y, ~, info] = EFGP(x, meas, sigma^2, ker, );
  fprintf('CPU times (s):'); fprintf('\t%.3g',info.cputime); fprintf('\n');
  fprintf('y.mean: rms resid of lin sys   %.3g\n', rms(y.mean-y.meanbook))
  fprintf('        rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
  % estim ability to average away noise via # pts in the rough kernel support...
  fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))

  % show pics
  if dim==1, figure; plot(x,meas,'.'); hold on; plot(x,y.mean,'-');
  elseif dim==2, figure;
    subplot(1,2,1); scatter(x(1,:),x(2,:),[],meas,'filled');
    caxis([-1 1]); axis equal tight
    subplot(1,2,2); scatter(x(1,:),x(2,:),[],y.mean,'filled');
    caxis([-1 1]); axis equal tight
  end
  title(sprintf('naive\\_gp test %dd',dim)); drawnow;
    
end             % ..........
