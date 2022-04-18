function [y, ytrg, info] = naive_gp(x, meas, sigmasq, ker, xtrg, opts)
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
%         getcovar - if true, compute conditional covariances [*** to implement]
%
% Outputs:
%  y - struct with fields of regression results corresp to given data points x:
%     mean - posterior mean vector, N*1
%     alpha - N*1 vector of weights for reconstruction of mean
%     meanbook - posterior mean vector, N*1, computed by "book formula" K*alpha
%     covar - posterior covariance, N*N. [optional]
%  ytrg - [optional; otherwise empty] struct of regression at new targets xtrg:
%     mean - posterior mean vector, n*1, computed by book formula
%     covar - posterior covariance, n*n. [optional]
%  info - diagnostic struct containing fields:
%     cputime - list of times in seconds for fill, solve + y eval, trg eval, etc.
%
% If called without arguments, does a self-test
if nargin==0, test_naive_gp; return; end
do_trg = (nargin>=5 && ~isempty(xtrg));
[dim,N] = size(x);
if numel(meas)~=N, error('sizes of meas and x must match!'); end
if N>1e4, warning('N getting too big for naive method!'); end

tic; K = densekermat(ker.k,x);
info.cputime(1) = toc;

tic;
meas = meas(:);
alpha = (K + sigmasq*eye(N)) \ meas;       % dense solve, probably
y.mean = meas - sigmasq*alpha;             % magic stable cheap formula
y.meanbook = K*alpha;                      % the book formula, slow, unstable
info.cputime(2) = toc;

ytrg = [];
if do_trg                              % *** untested
  [~,n] = size(xtrg);        % # targs
  tic;
  B = densekermat(ker.k,xtrg,x);   %  n-by-N
  ytrg.mean = B * alpha;    % do the GP sum, naively
  info.cputime(3) = toc;
end


%%%%%%%%%%
function test_naive_gp   % basic tests for now
N = 3e3;        % problem size
l = 0.1;        % SE kernel scale
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified

for dim = 1:3   % ..........
  fprintf('\ntest naive_gp, sigma=%.3g, dim=%d...\n',sigma,dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  ker = SE_ker(dim,l);

  ntrgs = 20;
  xtrg = equispaced_grid(dim, ntrgs);
  [y, ytrg, info] = naive_gp(x, meas, sigma^2, ker, xtrg);
  fprintf('CPU times (s):'); fprintf('\t%.3g',info.cputime); fprintf('\n');
  fprintf('y.mean: rms resid of lin sys   %.3g\n', rms(y.mean-y.meanbook))
  fprintf('        rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
  % estim ability to average away noise via # pts in the rough kernel support...
  fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))

  % show pics
  if dim==1, figure; plot(x,meas,'.'); hold on; plot(x,y.mean,'-');
      plot(xtrg,ytrg.mean,'-r');
  elseif dim==2, figure;
    subplot(1,2,1); scatter(x(1,:),x(2,:),[],meas,'filled');
    caxis([-1 1]); axis equal tight
    subplot(1,2,2); scatter(x(1,:),x(2,:),[],y.mean,'filled');
    caxis([-1 1]); axis equal tight
  end
  title(sprintf('naive\\_gp test %dd',dim)); drawnow;
    
end             % ..........
