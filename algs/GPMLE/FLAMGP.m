function [y, info] = FLAMGP(x, meas, sigmasq, ker, xtrg, opts)
% EFGP   GP regression via equispaced Fourier iterative method, in dim=1,2 or 3
%
% [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrg, opts)
%  performs Gaussian process regression using an arbitrary smooth(ish)
%  translation-invariant isotropic prior kernel ker (via a fast direct 
%  solver used to invert the covariance matrix),
%  conditioned on the y-values meas at the set of points x, in
%  spatial dimension 1,2 or 3. FLAM library is a prerequisite -- https://github.com/klho/FLAM
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
%         occ - max number of points per box, (default 64)
%         p - number of proxy points;
%         v - verbose;
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

% Notes: 1) target means currently not implemented
if nargin==0, test_FLAMGP; return; end

[dim,N] = size(x);
if nargin<5, xtrg = []; end
do_trg = ~isempty(xtrg);
if nargin<6, opts = []; end
if ~isfield(opts,'tol'), opts.tol = 1e-6; end     % default
if ~isfield(opts,'occ')
    if(dim == 1)
        opts.occ = 20;
    elseif(dim == 2)
        opts.occ = 64;
    end
end
if ~isfield(opts,'p'), opts.p = 16; end

if numel(meas)~=N, error('sizes of meas and x must match!'); end



verb = 1;

if(dim == 1)
    proxy = linspace(1.5,2.5,opts.p); proxy = [-proxy proxy];  % proxy points
elseif(dim == 2)
    theta_proxy = (1:opts.p)*2*pi/opts.p;
    proxy_ = [cos(theta_proxy); sin(theta_proxy)];
    proxy = [];
    for r = linspace(1.5,2.5,opts.p)
     proxy = [proxy r*proxy_];
    end
end

clear theta_proxy proxy_;

Afun = @(i,j)Afunflam(i,j,x,ker,sigmasq);
pxyfun = @(x,slf,nbr,l,ctr)pxyfunflam(x,slf,nbr,l,ctr,proxy,ker);


% verbose mode?
if (isfield(opts,'v') && (opts.v == true))
    opts_use = struct('symm','p','verb',verb);
else
    opts_use = struct('symm', 'p');
end

tt1 = tic;
F = rskelf(Afun,x,opts.occ,opts.tol,pxyfun,opts_use);
alpha = rskelf_sv(F,meas);
y.mean = meas - alpha*sigmasq;
info.cpu_time = toc(tt1);


    



%%%%%%%%%%
function test_FLAMGP   % basic tests for now, duplicates naive_gp *** to unify
N = 3e3;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-8;

for dim = 1:2   % ..........
  fprintf('\ntest EFGP, sigma=%.3g, tol=%.3g, dim=%d...\n',sigma,opts.tol,dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  ker = SE_ker(dim,l);
  [y, ~] = FLAMGP(x, meas, sigma^2, ker, [], opts);
  % run o(n^3) naive gp regression
  [ytrue, ytrg, ~] = naive_gp(x, meas, sigma^2, ker, [], opts);
  %fprintf('%d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
  %fprintf('CPU times (s):'); fprintf('\t%.3g',info.cputime); fprintf('\n');
  fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
  % estim ability to average away noise via # pts in the rough kernel support...
  fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))
  % make sure we're computing gp regression accurately
  fprintf('        rms flamgp vs naive      %.3g\n', rms(y.mean-ytrue.mean))

  % show pics
  if dim==1, figure; plot(x,meas,'.'); hold on; plot(x,y.mean,'-');
  elseif dim==2, figure;
    subplot(1,2,1); scatter(x(1,:),x(2,:),[],meas,'filled');
    caxis([-1 1]); axis equal tight
    subplot(1,2,2); scatter(x(1,:),x(2,:),[],y.mean,'filled');
    caxis([-1 1]); axis equal tight
  end
  title(sprintf('FLAMGP test %dd',dim)); drawnow;

end % ..........
