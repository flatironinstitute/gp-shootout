function [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrg, opts)
% FLAMGP  GP regression via recursive skeletonization, in dim 1 or 2.
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
%         tol - desired tolerance, default: 1e-6
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
%     cpu_time - run time in seconds
%     RAM - memory usage in bytes
%     proxy - list of proxy pts used (for unit scale box)
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
if ~isfield(opts,'p'), opts.p = ceil(log(opts.tol)/log(sqrt(2.0)/3.0)/2)  ; end

fprintf('p = %d\n',opts.p);

if numel(meas)~=N, error('sizes of meas and x must match!'); end

verb = 1;

if(dim == 1)
    proxy = linspace(1.5,2.5,opts.p); proxy = [-proxy proxy];  % proxy points
elseif(dim == 2)
    theta_proxy = (1:opts.p)*2*pi/opts.p;
    proxy_ = [cos(theta_proxy); sin(theta_proxy)];
    proxy = [];
    for r = linspace(1.5,2.5,opts.p)        % mysterious shells
     proxy = [proxy r*proxy_];
    end
else
  error('dim>2 not implemented!');
end

clear theta_proxy proxy_;

Afun = @(i,j) Afunflam(i,j,x,ker,sigmasq);
pxyfun = @(x,slf,nbr,l,ctr) pxyfunflam(x,slf,nbr,l,ctr,proxy,ker);


% verbose mode?
if (isfield(opts,'v') && (opts.v == true))
    opts_use = struct('symm','p','verb',opts.v);
else
    opts_use = struct('symm', 'p');
end

tt1 = tic;
F = rskelf(Afun,x,opts.occ,opts.tol,pxyfun,opts_use);
alpha = rskelf_sv(F,meas);         % soln vec via solve the lin sys
y.mean = meas - alpha*sigmasq;     % magic formula for means at data pts
info.cpu_time = toc(tt1);
info.proxy = proxy;
w = whos('F');
info.RAM = w.bytes;
dummy = [];

% posterior mean at targets
ytrg = [];
if do_trg
   tic;
%   [ndim,ntrg] = size(xtrg);        % # targs
%   ntot = N + ntrg;
%   xflam_targ = zeros(ndim,ntot);
%   xflam_targ(:,1:N) = x;
%   xflam_targ(:,N+1:end) = xtrg;
  %Afun_targ = @(i,j) Afunflam_targ(i,j,xflam_targ,ker,N);
  Afun_targ = @(i,j) Afunflam_targ(i,j,xtrg,x,ker);
  %pxyfun_targ = @(x,slf,nbr,l,ctr) pxyfunflam_targ(xflam_targ,slf,nbr,l,ctr,proxy,ker,N);
  pxyfun_targ = @(rc,xtrg,x,slf,nbr,l,ctr) pxyfunflam_targ(rc,xtrg,x,slf,nbr,l,ctr,proxy,ker);
  
  if (isfield(opts,'v') && (opts.v == true))
    opts_use = struct('symm','n','verb',verb);
  else
    opts_use = struct('symm', 'n');
  end
  warning('off');
  F_targ = rskel(Afun_targ,xtrg,x,opts.occ,opts.tol,pxyfun_targ,opts_use);
  warning('on');
%   alpha_use = zeros(ntot,1);
%   alpha_use(1:N) = alpha;
  ytrg.mean = rskel_mv(F_targ,alpha);
  
  %B = densekermat(ker.k,xtrg,x);   %  n-by-N
  %ytrg.mean = ytrg_tmp(N+1:end);    % do the GP sum, naively
  info.cpu_time(3) = toc;
end




%%%%%%%%%%
function test_FLAMGP   % basic tests for now, duplicates naive_gp *** to unify
N = 3e3;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-10;



for dim = 2:2   % ..........
  fprintf('\ntest EFGP, sigma=%.3g, tol=%.3g, dim=%d...\n',sigma,opts.tol,dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  xtrg = rand(dim,N);
  ker = SE_ker(dim,l);
  [y, ytrg, info] = FLAMGP(x, meas, sigma^2, ker, xtrg, opts);
  % run O(n^3) naive gp regression
  [ytrue, ytrg_true, ~] = naive_gp(x, meas, sigma^2, ker, xtrg, opts);
  fprintf('%d proxies \t %.g GB RAM\n',numel(info.proxy),info.RAM/1e9)
  fprintf('CPU time (s):'); fprintf('\t%.3g',info.cpu_time); fprintf('\n');
  fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
  % estim ability to average away noise via # pts in the rough kernel support...
  fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))
  % make sure we're computing gp regression accurately
  fprintf('        rms flamgp vs naive      %.3g\n', rms(y.mean-ytrue.mean))
  fprintf('        rms flamgp vs naive at targets    %.3g\n', rms(ytrg.mean-ytrg_true.mean))

  % show pics
  if dim==1, figure; plot(x,meas,'.'); hold on; plot(x,y.mean,'-');
  elseif dim==2, figure;
    subplot(1,2,1); scatter(x(1,:),x(2,:),[],meas,'filled');
    caxis([-1 1]); axis equal tight
    subplot(1,2,2); scatter(x(1,:),x(2,:),[],y.mean,'filled');
    caxis([-1 1]); axis equal tight
  end
  title(sprintf('FLAMGP test %dd',dim)); drawnow;
end