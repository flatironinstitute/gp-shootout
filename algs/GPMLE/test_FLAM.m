N = 1e5;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale
sigma = 0.5;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 10.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified

opts.v = true;

opts.tol = 1e-10;

for dim = 2:2   % ..........
  fprintf('\ntest EFGP, sigma=%.3g, tol=%.3g, dim=%d...\n',sigma,opts.tol,dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  xtrg = rand(dim,ceil(N/100));
  
  ntrgs_per_dim = 100;
  xtrgs = equispaced_grid(dim, ntrgs_per_dim);
  xmax = max(x, [], 'all');
  xmin = min(x, [], 'all');
  xtrgs = xtrgs .* (xmax - xmin) + xmin;

  ker = SE_ker(dim,l);
  nu = 0.5;
%  ker = Matern_ker(dim, nu, l);
  tic, [y, ytrg, info] = FLAMGP(x, meas, sigma^2, ker, xtrg, opts); toc

  return
  
  
  opts.tol = 1e-12;
  [ytrue, ytrg_true, ~] = FLAMGP(x, meas, sigma^2, ker, xtrg, opts);

  
%   % run O(n^3) naive gp regression
%   [ytrue, ytrg_true, ~] = naive_gp(x, meas, sigma^2, ker, xtrg, opts);
  fprintf('%d proxies \t %.g GB RAM\n',numel(info.proxy),info.RAM/1e9)
  fprintf('CPU time (s):\n'); disp(info.cpu_time)
  
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
warning('on')
