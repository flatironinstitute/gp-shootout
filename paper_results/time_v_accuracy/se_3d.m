% accuracy vs time for for fixed N
clear opts;

dim = 3;

% set directory for saving results and loading data
dir = "~/gp-shootout/results/philip/time_v_accuracy/data";

% sigma used to generate data and to be used for regression
load(fullfile(dir, 'sigmatrue.mat'));

% load data
load(fullfile(dir, 'x_3d_1e5.mat'));
load(fullfile(dir, 'meas_3d_1e5.mat'));

% targets
ntrgs_per_dim = 30;
xtrgs = equispaced_grid(dim, ntrgs_per_dim);
xmax = max(x, [], 'all');
xmin = min(x, [], 'all');
xtrgs = xtrgs .* (xmax - xmin) + xmin;
opts.only_trgs = 1;

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
ker = SE_ker(dim,l);



%x = x(:,1:100);
%meas = meas(1:100);
%[y1, ytrg1, info1] = naive_gp(x, meas, sigmasq, ker, xtrgs, []);

% % Accurate solution via efgp
% opts.tol = 1e-10;
% [y_true0, ytrg_true0, info0] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
% fprintf('\n');
% 
% opts.tol = 1e-12;
% [y_true, ytrg_true, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
% fprintf('\n');
% save(fullfile(dir, 'se_3d_true.mat'), 'ytrg_true');
% save(fullfile(dir, 'se_3d_true0.mat'), 'ytrg_true0');

load(fullfile(dir, 'se_3d_true.mat'));
load(fullfile(dir, 'se_3d_true0.mat'));

fprintf('max dd: %g\n', max(abs(ytrg_true0.mean - ytrg_true.mean)));
fprintf('rms dd: %g\n', rms(ytrg_true0.mean - ytrg_true.mean));


% EFGP
nns = 5;
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

efgp_3d.ts = ts;
efgp_3d.rms_errs = rms_errs;
efgp_3d.linf_errs = linf_errs;
save(fullfile(dir, 'efgp_3d.mat'), 'efgp_3d');




%RLCM
nns = 3;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.rank = 50 * i;
    [y, ytrg, info] = RLCM(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

rlcm_3d.ts = ts;
rlcm_3d.rms_errs = rms_errs;
rlcm_3d.linf_errs = linf_errs;
save(fullfile(dir, 'rlcm_3d.mat'), 'rlcm_3d');




% SKI
nns = 5;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    fprintf('%g\n', i);
    opts.grid_size = 10 * i;
    [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrgs, opts);

    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end
ski_3d.ts = ts;
ski_3d.rms_errs = rms_errs;
ski_3d.linf_errs = linf_errs;
save(fullfile(dir, 'ski_3d.mat'), 'ski_3d');





% FLAM
opts.v = true;
    
nns = 2;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.tol = 1e-3 * 10^(-i);
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_3d.ts = ts;
flam_3d.rms_errs = rms_errs;
flam_3d.linf_errs = linf_errs;
save(fullfile(dir, 'flam_3d.mat'), 'flam_3d');




% FLAM no proxy
opts.v = true;
opts.no_proxy = true;

nns = 2;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.tol = 1e-3 * 10^(-i);
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_nopxy_3d.ts = ts;
flam_nopxy_3d.rms_errs = rms_errs;
flam_nopxy_3d.linf_errs = linf_errs;
save(fullfile(dir, 'flam_nopxy_3d.mat'), 'flam_nopxy_3d');




% plotting
load(fullfile(dir, 'efgp_3d.mat'));
load(fullfile(dir, 'ski_3d.mat'));
load(fullfile(dir, 'flam_3d.mat'));
load(fullfile(dir, 'flam_nopxy_3d.mat'));
load(fullfile(dir, 'rlcm_3d.mat'));
hold on;
plot(log10(efgp_3d.rms_errs), log10(efgp_3d.ts), '-o');
plot(log10(rlcm_3d.rms_errs), log10(rlcm_3d.ts), '-o');
plot(log10(ski_3d.rms_errs), log10(ski_3d.ts), '-o');
plot(log10(flam_3d.rms_errs), log10(flam_3d.ts), '-o');
plot(log10(flam_nopxy_3d.rms_errs), log10(flam_3d.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;


% print for paper
fprintf('EFGP\n')
print_tikz(log10(efgp_3d.rms_errs), log10(efgp_3d.ts))
fprintf('\nSKI\n')
print_tikz(log10(ski_3d.rms_errs), log10(ski_3d.ts))
fprintf('\nFLAM\n')
print_tikz(log10(flam_3d.rms_errs), log10(flam_3d.ts))
fprintf('\nFLAM no pxy\n')
print_tikz(log10(flam_nopxy_3d.rms_errs), log10(flam_nopxy_3d.ts))
fprintf('\nRLCM\n')
print_tikz(log10(rlcm_3d.rms_errs), log10(rlcm_3d.ts))
