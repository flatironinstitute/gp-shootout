% this script compares the accuracy of several algorithms for 
% Gaussian process regression against the naive, O(N^3) algorithm.
% specifcally, this script compares the posterior mean at the inputted data
% points using several fast methods and compares that to the same values
% obtained using a slow, accurate method.
%
% It is thus part of the accuracy experiments rather than a self-test for
% any particular method (AHB), hence not in test/
%
% Note: tol is now different for each method, more useful.

N = 5000;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale
sigma = 0.03;    % used to regress (1=easy; 0.01 = ill cond kills EFGP,SKI)
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified


for dim = 1:3

  fprintf('\ndim=%g...\n', dim)
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
  ker = SE_ker(dim,l);        % SE needed for RLCM for now
  
  % run O(n^3) naive gp regression (backward stable)
  [ytrue, ytrg, ~] = naive_gp(x, meas, sigma^2, ker, [], []);
 
  % EFGP
  opts.tol = 1e-6;
  [y1, ~, info1] = EFGP(x, meas, sigma^2, ker, [], opts);
  fprintf('EFGP   rms vs naive %.3g,\t time: %.3g s\n', rms(y1.mean-ytrue.mean), info1.cpu_time.total);
  
  % FLAMGP (currently throws error in 3d)
  if dim < 3
    opts = [];
    opts.tol = 1e-12;      % crucial: too big and get chol err!
    [y3, ~, info3] = FLAMGP(x, meas, sigma^2, ker, [], opts);
    fprintf('FLAMGP rms vs naive %.3g,\t time: %.3g s\n', rms(y3.mean-ytrue.mean), info3.cpu_time.total);
  end

  % SKI w/ default gridsize
  opts = [];             % start afresh. There's no tol param for SKI
  opts.grid_size = 100;   % 1000?
  [y2, ~, info2] = SKI(x, meas, sigma^2, ker, [], opts);
  fprintf('SKI    rms vs naive %.3g,\t time: %.3g s\n', rms(y2.mean-ytrue.mean), info2.cpu_time.total);

  % RLCM - needs user to have compiled, via, eg algs/RCLM/buildit.sh...
  opts = [];       % there's no tol param for RLCM
  opts.rank = 200;      % some kinda convergence param (large needed in 3D?)
  [y4, ~, info4] = RLCM(x, meas, sigma^2, ker, [], opts);
  fprintf('RLCM   rms vs naive %.3g,\t time: %.3g s\n', rms(y4.mean-ytrue.mean), info4.cpu_time.total);

end
