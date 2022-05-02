% simply try to run all meths on a larger example. AHB 5/2/22
clear; verb = 0;

N = 1e5;        % problem size  (>=1e5 breaks SKI)
l = 0.1;        % kernel scale
sigma = 0.1;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-6;   % not used by all meths
dim = 2;
unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
wavevec = freqdata*unitvec;    % col vec
f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
[x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);

ker = SE_ker(dim,l);
%ker = Matern_ker(dim,3/2,l);   % harder for EFGP; not yet in RLCM

fprintf('\ntest N=%d, sigma=%.3g, tol=%.3g, dim=%d...\n',N,sigma,opts.tol,dim)

% *** to do: compare soln vectors!
for meth=[1 2 4]
  switch meth
    case 1, disp('EFGP...')
      [y, ~, info] = EFGP(x, meas, sigma^2, ker, [], opts);     % regress
      fprintf('EFGP %d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
    case 2, disp('FLAMGP...')
      [y, ~, info] = FLAMGP(x, meas, sigma^2, ker, [], opts);     % 1 min for N=1e6
      fprintf('%d proxies \t %.g GB RAM, times:\n',numel(info.proxy),info.RAM/1e9)
    case 3, disp('SKI...')       % *** has problems even passing back small N
      [y, ~, info] = SKI(x, meas, sigma^2, ker, [], opts);
    case 4, disp('RLCM...')
      [y, ~, info] = RLCM(x, meas, sigma^2, ker, [], opts);
  end
  fprintf('CPU time (s):\n'); disp(info.cpu_time);
  fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
  % estim ability to average away noise via # pts in the rough kernel support...
  fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))
end
  
if verb, figure; siz = 1.0; Np = min(N,1e5);  % max pts to plot
  subplot(1,2,1); scatter(x(1,1:Np),x(2,1:Np),siz*ones(Np,1),meas(1:Np),'filled');
  caxis([-1 1]); axis equal tight
  subplot(1,2,2); scatter(x(1,1:Np),x(2,1:Np),siz*ones(Np,1),y.mean(1:Np),'filled');
  caxis([-1 1]); axis equal tight
  title(sprintf('expt: %dd N=%d',dim,N));
end
