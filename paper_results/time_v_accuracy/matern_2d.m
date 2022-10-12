% accuracy vs time for for fixed N
clear opts;

dim = 2;

% set directory for saving results and loading data
dir = "~/gp-shootout/results/philip/time_v_accuracy/data";

% sigma used to generate data and to be used for regression
load(fullfile(dir, 'sigmatrue.mat'));

% load data
load(fullfile(dir, 'x_2d_1e5.mat'));
load(fullfile(dir, 'meas_2d_1e5.mat'));

% targets
ntrgs = 100;
xtrgs = equispaced_grid(dim, ntrgs);
opts.only_trgs = 1;

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
var = 1;
nu = 0.5;
ker = Matern_ker(dim, nu, l, var);

% get accurate solution
% opts.v = true;
% opts.tol = 1e-8;
% [y0, ytrg_true0, info0] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
% opts.tol = 1e-10;
% [y, ytrg_true, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
% save(fullfile(dir, 'matern_2d_true.mat'), 'ytrg_true');
% save(fullfile(dir, 'matern_2d_true0.mat'), 'ytrg_true0');

load(fullfile(dir, 'matern_2d_true.mat'));
load(fullfile(dir, 'matern_2d_true0.mat'));
fprintf('max dd: %g\n', max(abs(ytrg_true.mean - ytrg_true0.mean)));




% EFGP
nns = 3;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.tol = 1e-1 * 10^(-i);
    [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
    
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

efgp_2d_matern.ts = ts;
efgp_2d_matern.rms_errs = rms_errs;
efgp_2d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'efgp_2d_matern.mat'), 'efgp_2d_matern');


% SKI
nns = 3;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    fprintf('i: %d\n', i);
    opts.grid_size = 10^i;
    [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrgs, opts);

    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

ski_2d_matern.ts = ts;
ski_2d_matern.rms_errs = rms_errs;
ski_2d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'ski_2d_matern.mat'), 'ski_2d_matern');



% FLAM
nns = 5;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
opts.v = true;
for i=1:nns
    disp(i);
    opts.tol = 10^(-4 - i); 
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_2d_matern.ts = ts;
flam_2d_matern.rms_errs = rms_errs;
flam_2d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'flam_2d_matern.mat'), 'flam_2d_matern');




%RLCM
nns = 3;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.rank = 100 * i;
    [y, ytrg, info] = RLCM(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

rlcm_2d_matern.ts = ts;
rlcm_2d_matern.rms_errs = rms_errs;
rlcm_2d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'rlcm_2d_matern.mat'), 'rlcm_2d_matern');





% plotting
load(fullfile(dir, "efgp_2d_matern.mat"));
load(fullfile(dir, "ski_2d_matern.mat"));
load(fullfile(dir, "flam_2d_matern.mat"));
load(fullfile(dir, "rlcm_2d_matern.mat"));
hold on;
plot(log10(efgp_2d_matern.rms_errs), log10(efgp_2d_matern.ts), '-o');
plot(log10(ski_2d_matern.rms_errs), log10(ski_2d_matern.ts), '-o');
plot(log10(flam_2d_matern.rms_errs), log10(flam_2d_matern.ts), '-o');
plot(log10(rlcm_2d_matern.rms_errs), log10(rlcm_2d_matern.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;




% print for paper
fprintf('EFGP\n')
print_tikz(log10(efgp_2d_matern.rms_errs), log10(efgp_2d_matern.ts))
fprintf('\nSKI\n')
print_tikz(log10(ski_2d_matern.rms_errs), log10(ski_2d_matern.ts))
fprintf('\nFLAM\n')
print_tikz(log10(flam_2d_matern.rms_errs), log10(flam_2d_matern.ts))
fprintf('\nRLCM\n')
print_tikz(log10(rlcm_2d_matern.rms_errs), log10(rlcm_2d_matern.ts))