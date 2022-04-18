function [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrg, opts)
% EFGP   GP regression via equispaced Fourier iterative method, in dim=1,2 or 3
%
% [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrg, opts)
%  performs Gaussian process regression using an arbitrary smooth(ish)
%  translation-invariant isotropic prior kernel ker (via its Fourier transform
%  ker.khat), conditioned on the y-values meas at the set of points x, in
%  spatial dimension 1,2 or 3. The method is efficient and accurate only when
%  khat decays rapidly in Fourier space. FINUFFT library is a prerequisite.
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
%  ytrg - [optional; otherwise empty] struct of regression at new targets xtrg:
%     mean - posterior mean vector, n*1
%  info - diagnostic struct containing fields:
%     xis - Fourier xi nodes use
%     beta - m*1 vector of weight-space (Fourier basis) weights
%     cputime - list of times in seconds for various steps
%     iter - # iterations needed
%
% If called without arguments, does a self-test.

% Notes: 1) this code is a wrapper to separate dimension functions
%   2) I changed the ker format from Philip
% Todo:
% *** figure how to switch off x as self-targs, since may make too slow?
if nargin==0, test_EFGP; return; end
if nargin<5, xtrg = []; end
do_trg = ~isempty(xtrg);
if nargin<6, opts = []; end
if ~isfield(opts,'tol'), opts.tol = 1e-6; end     % default

[dim,N] = size(x);
if numel(meas)~=N, error('sizes of meas and x must match!'); end
n = size(xtrg,2);   % # new targets

xsol = [x, xtrg]';  % hack for now which adds meas pts to target list
                    % and transpose to Philip n*d shape                    
if dim==1
  [info.beta, info.xis, yhat, info.iter, info.cputime] = function_space1d(x', meas, sigmasq, ker, opts.tol, xsol);
elseif dim==2
  [info.beta, info.xis, yhat, info.iter, info.cputime] = function_space2d(x', meas, sigmasq, ker, opts.tol, xsol); 
elseif dim==3
  [info.beta, info.xis, yhat, info.iter, info.cputime] = function_space3d(x', meas, sigmasq, ker, opts.tol, xsol); 
else
  error('dim must be 1,2, or 3!');
end

y.mean = yhat(1:N);   % hack for now to split out posterior means into two types
ytrg.mean = yhat(N+1:end);


%%%%%%%%%%
function test_EFGP   % basic tests for now, duplicates naive_gp *** to unify
N = 3e3;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-8;

for dim = 1:3   % ..........
  fprintf('\ntest EFGP, sigma=%.3g, tol=%.3g, dim=%d...\n',sigma,opts.tol,dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  ker = SE_ker(dim,l);
  [y, ~, info] = EFGP(x, meas, sigma^2, ker, [], opts);
  % run o(n^3) naive gp regression
  [ytrue, ytrg, ~] = naive_gp(x, meas, sigma^2, ker, [], opts);
  fprintf('%d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
  fprintf('CPU times (s):'); fprintf('\t%.3g',info.cputime); fprintf('\n');
  fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
  % estim ability to average away noise via # pts in the rough kernel support...
  fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))
  % make sure we're computing gp regression accurately
  fprintf('        rms efgp vs naive      %.3g\n', rms(y.mean-ytrue.mean))

  % show pics
  figure;
  if dim==1, plot(x,meas,'.'); hold on; plot(x,y.mean,'-');
  elseif dim==2
    subplot(1,2,1); scatter(x(1,:),x(2,:),[],meas,'filled');
    caxis([-1 1]); axis equal tight
    subplot(1,2,2); scatter(x(1,:),x(2,:),[],y.mean,'filled');
    caxis([-1 1]); axis equal tight
  elseif dim==3
    subplot(1,2,1); scatter3(x(1,:),x(2,:),x(3,:),[],meas,'filled');
    caxis([-1 1]); axis equal tight
    subplot(1,2,2); scatter3(x(1,:),x(2,:),x(3,:),[],y.mean,'filled');
    caxis([-1 1]); axis equal tight
  end
  title(sprintf('EFGP test %dd',dim)); drawnow;

end             % ..........
