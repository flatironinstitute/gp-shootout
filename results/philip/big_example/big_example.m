rng(1);

% set directory for loading data and saving results
dir = "~/gp-shootout/results/philip/big_example/data";
load(fullfile(dir, 'sigmatrue.mat'));

% set tolerances and make sure we only evaluate posterior mean at trgs
opts2d.tol = 1e-5;
opts2d.only_trgs = true;
opts2d_true.tol = 1e-6;
opts2d_true.only_trgs = true;

% load data
load(fullfile(dir, 'x_2d.mat'));
load(fullfile(dir, 'meas_2d.mat'));

% subsample data
N = 3e8;
meas = meas(1:N);
x = x(:,1:N);

% gp regression
dim = 2;
l = 0.3;
nu = 1.5;
var = 1.0;
ker = Matern_ker(dim, nu, l, var);
%ker = SE_ker(dim, l, var);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);

sigmasq = sigmatrue^2;
[y, ytrgs, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts2d);
[y_true, ytrgs_true, info_true] = EFGP(x, meas, sigmasq, ker, xtrgs, opts2d_true);

err_rms = rms(ytrgs.mean-ytrgs_true.mean);
err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));

m = (numel(info.xis) -  1) / 2;
fprintf("$ 10^{%d}$ & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
     log10(N), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);

filename = sprintf('mat_%gd_info_1e%g.mat', dim, log10(N));
save(fullfile(dir, filename), 'info')
filename = sprintf('mat_%gd_rms_err_1e%g.mat', dim, log10(N));
save(fullfile(dir, filename), 'err_rms')


