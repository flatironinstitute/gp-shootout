% this script compares the accuracy of several algorithms for 
% Gaussian process regression against the naive, O(N^3) algorithm.
% specifcally, this script compares the posterior mean at the inputted data
% points using several fast methods and compares that to the same values
% obtained using a slow, accurate method. 

N = 30;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified



for dim = 1:3
  fprintf('dim=%g...\n', dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  ker = SE_ker(dim,l);
  % run o(n^3) naive gp regression
  [ytrue, ytrg, ~] = naive_gp(x, meas, sigma^2, ker, [], []);
 
  % EFGP
  opts.tol = 1e-8;
  [y1, ~, info] = EFGP(x, meas, sigma^2, ker, [], opts);
  fprintf('EFGP   rms vs naive %.3g, time: %.3g\n', rms(y1.mean-ytrue.mean), info.cpu_time(end));

  % SKI w/ default gridsize.    *** does it need funky x=0,1 pts present?
  %%%opts.grid_size = 1000;
  [y2, ~, info] = SKI(x, meas, sigma^2, ker, [], opts);
  fprintf('SKI    rms vs naive %.3g, time: %.3g\n', rms(y2.mean-ytrue.mean), info.cpu_time(end));

  % FLAMGP - currently throws error in 3d
  if dim < 3 
      [y3, ~, info] = FLAMGP(x, meas, sigma^2, ker, [], opts);
      fprintf('FLAMGP rms vs naive %.3g, time: %.3g\n', rms(y3.mean-ytrue.mean), info.cpu_time(end));
  end
end
