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
ntrgs_per_dim = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_dim);
xmax = max(x, [], 'all');
xmin = min(x, [], 'all');
xtrgs = xtrgs .* (xmax - xmin) + xmin;
opts.only_trgs = 1;

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
dim = 2;
ker = SE_ker(dim,l);

opts.tol = 1e-8;
[y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
opts.tol = 1e-10;
[y2, ytrg2, info2] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));



% get accurate solution
opts.tol = 1e-12;
[y_true0, ytrg_true0, info0] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
fprintf('\n');

opts.tol = 1e-14;
[y_true, ytrg_true, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
fprintf('\n');
%%%scatter(x, meas); hold on; plot(xtrgs, ytrg_true.mean); hold off
%%%[y, ytrg_true, ~] = naive_gp(x, meas, sigmasq, ker, xtrgs, []);
fprintf('max dd: %g\n', max(abs(ytrg_true0.mean - ytrg_true.mean)));
fprintf('rms dd: %g\n', rms(ytrg_true0.mean - ytrg_true.mean));




%RLCM
nns = 4;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.rank = 10 * i;
    [y, ytrg, info] = RLCM(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

rlcm_2d.ts = ts;
rlcm_2d.rms_errs = rms_errs;
rlcm_2d.linf_errs = linf_errs;
save(fullfile(dir, 'rlcm_2d.mat'), 'rlcm_2d');



% EFGP
nns = 6;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.tol = 1e-2 * 10^(-i);
    [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

efgp_2d.ts = ts;
efgp_2d.rms_errs = rms_errs;
efgp_2d.linf_errs = linf_errs;
save(fullfile(dir, 'efgp_2d.mat'), 'efgp_2d');

%disp(efgp_2d.rms_errs);
%plot(log10(efgp_2d.rms_errs), log10(efgp_2d.ts), '-o');
%set ( gca, 'xdir', 'reverse' );

% SKI
nns = 2;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.grid_size = 10^i;
    [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrgs, opts);

    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

ski_2d.ts = ts;
ski_2d.rms_errs = rms_errs;
ski_2d.linf_errs = linf_errs;
save(fullfile(dir, 'ski_2d.mat'), 'ski_2d');



% FLAM
nns = 4;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.tol = 1e-6 * 10^(-i);
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_2d.ts = ts;
flam_2d.rms_errs = rms_errs;
flam_2d.linf_errs = linf_errs;
save(fullfile(dir, 'flam_2d.mat'), 'flam_2d');





% plotting
load(fullfile(dir, 'efgp_2d.mat'));
load(fullfile(dir, 'ski_2d.mat'));
load(fullfile(dir, 'flam_2d.mat'));
load(fullfile(dir, 'rlcm_2d.mat'));
hold on;
plot(log10(efgp_2d.rms_errs), log10(efgp_2d.ts), '-o');
plot(log10(ski_2d.rms_errs), log10(ski_2d.ts), '-o');
plot(log10(flam_2d.rms_errs), log10(flam_2d.ts), '-o');
plot(log10(rlcm_2d.rms_errs), log10(rlcm_2d.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;


% print for paper
fprintf('EFGP\n')
print_tikz(log10(efgp_2d.rms_errs), log10(efgp_2d.ts))
fprintf('\nSKI\n')
print_tikz(log10(ski_2d.rms_errs), log10(ski_2d.ts))
fprintf('\nFLAM\n')
print_tikz(log10(flam_2d.rms_errs), log10(flam_2d.ts))
fprintf('\nRLCM\n')
print_tikz(log10(rlcm_2d.rms_errs), log10(rlcm_2d.ts))
