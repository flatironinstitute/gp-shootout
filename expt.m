% driver script for an experiment
% (copied from self-test of EFGP)
clear; verb = 1;

N = 1e6;        % problem size
l = 0.1;        % SE kernel scale
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-8;
dim = 2;
unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
wavevec = freqdata*unitvec;    % col vec
f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
[x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
ker = SE_ker(dim,l);

fprintf('\ntest EFGP, N=%d, sigma=%.3g, tol=%.3g, dim=%d...\n',N,sigma,opts.tol,dim)
[y, ~, info] = EFGP(x, meas, sigma^2, ker, [], opts);     % regress

y.mean = real(y.mean);    % *** decide if complex-valued y is ok?
fprintf('%d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis),rms(info.beta))
fprintf('CPU times (s):'); fprintf('\t%.3g',info.cputime); fprintf('\n');
fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
% estim ability to average away noise via # pts in the rough kernel support...
fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))

if verb, figure; siz = 1.0;
  subplot(1,2,1); scatter(x(1,:),x(2,:),siz+0*meas,meas,'filled');
  caxis([-1 1]); axis equal tight
  subplot(1,2,2); scatter(x(1,:),x(2,:),siz+0*meas,y.mean,'filled');
  caxis([-1 1]); axis equal tight
  title(sprintf('expt: EFGP %dd N=%d',dim,N));
end