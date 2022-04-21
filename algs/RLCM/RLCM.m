function [y, ytrg, info] = RLCM(x, meas, sigmasq, ker, xtrg, opts)
% RLCM   GP regression via recursively Low-rank compressed matrices, Chen-Stein
%
% [y, ytrg, info] = RLCM(x, meas, sigmasq, ker, xtrg, opts)
%  performs Gaussian process regression using certain kernels.
%  *** finish.
%
%  Coordinates x and xtrg may be anywhere in R^2.
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
%         dense - use a (usually) slower dense linear solve instead of the
%         iterative solver 
%
% Outputs:
%  y - struct with fields of regression results corresp to given data points x:
%     mean - posterior mean vector, N*1
%  ytrg - [optional; otherwise empty] struct of regression at new targets xtrg:
%     mean - posterior mean vector, n*1
%  info - diagnostic struct containing fields:
%     cpu_time - struct with fields including (in seconds): total
%
% If called without arguments, does a self-test.

if nargin==0, test_RLCM; return; end
if nargin<5, xtrg = []; end
do_trg = ~isempty(xtrg);
if nargin<6, opts = []; end
if ~isfield(opts,'tol'), opts.tol = 1e-6; end     % default

[dim,N] = size(x);
if numel(meas)~=N, error('sizes of meas and x must match!'); end
n = size(xtrg,2);   % # new targets

tmpkey = randi([0 intmax('uint32')], 'uint32');
filetrain = sprintf('/tmp/RLCM_train_%x.tmp',tmpkey);
filetest = sprintf('/tmp/RLCM_test_%x.tmp',tmpkey);
fileymean = sprintf('/tmp/RLCM_ymean_%x.tmp',tmpkey);
fileytrgmean = sprintf('/tmp/RLCM_ytrgmean_%x.tmp',tmpkey);

fopen(

t0=tic;

***

info.cpu_time.total = toc(t0);




%%%%%%%%%%
function test_RLCM   % basic tests for now, duplicates naive_gp *** to unify
N = 3e3;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale rel to domain [0,1]^dim, ie hardness of prob
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-8;
L=0; shift = 0;
%L = 50.0; shift = 200;   % arbitary, tests correct centering and L-box rescale

for dim = 1:3   % ..........
  fprintf('\ntest RLCM, sigma=%.3g, tol=%.3g, dim=%d...\n',sigma,opts.tol,dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);    % x in [0,1]^dim
  x = L*x + (2*rand(dim,1)-1)*shift;                           % scale & shift
  ker = SE_ker(dim,L*l);             % note L also scales kernel length here
  [y, ~, info] = RLCM(x, meas, sigma^2, ker, [], opts);
  % run O(n^3) naive gp regression
  [ytrue, ytrg, ~] = naive_gp(x, meas, sigma^2, ker, [], opts);
  fprintf('%d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
  fprintf('CPU time (s):\t%.3g\n',info.cpu_time.total);
  fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
  % estim ability to average away noise via # pts in the rough kernel support...
  fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))
  % make sure we're computing gp regression accurately
  fprintf('        rms RLCM vs naive      %.3g\n', rms(y.mean-ytrue.mean))

  if 0       % show pics
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
    title(sprintf('RLCM test %dd',dim)); drawnow;
  end
  
end             % ..........
