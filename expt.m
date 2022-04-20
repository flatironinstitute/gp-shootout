% driver script for an experiment.
% (copied from self-test of EFGP)
clear; verb = 1;

N = 1e5;        % problem size
l = 0.1;        % kernel scale
sigma = 0.1;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-6;
dim = 2;
unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
wavevec = freqdata*unitvec;    % col vec
f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
[x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
ker = SE_ker(dim,l);
%ker = Matern_ker(dim,3/2,l);   % harder for EFGP

fprintf('\ntest N=%d, sigma=%.3g, tol=%.3g, dim=%d...\n',N,sigma,opts.tol,dim)
[y, ~, info] = EFGP(x, meas, sigma^2, ker, [], opts);     % regress
fprintf('EFGP %d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))

%[y, ~, info] = FLAMGP(x, meas, sigma^2, ker, [], opts);     % 1 min for N=1e6
%fprintf('%d proxies \t %.g GB RAM\n',numel(info.proxy),info.RAM/1e9)

fprintf('CPU times (s):'); fprintf('\t%.3g',info.cpu_time); fprintf('\n');
fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
% estim ability to average away noise via # pts in the rough kernel support...
fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))

if verb, figure; siz = 1.0; Np = min(N,1e5);  % max pts to plot
  subplot(1,2,1); scatter(x(1,1:Np),x(2,1:Np),siz*ones(Np,1),meas(1:Np),'filled');
  caxis([-1 1]); axis equal tight
  subplot(1,2,2); scatter(x(1,1:Np),x(2,1:Np),siz*ones(Np,1),y.mean(1:Np),'filled');
  caxis([-1 1]); axis equal tight
  title(sprintf('expt: %dd N=%d',dim,N));
end
