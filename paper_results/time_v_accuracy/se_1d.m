% accuracy vs time for for fixed N
clear opts;

dim = 1;

% set directory for saving results and loading data
dir = [fileparts(mfilename('fullpath')) '/data'];

% sigma used to generate data and to be used for regression
load(fullfile(dir, 'sigmatrue.mat'));

% load data
load(fullfile(dir, 'x_1d_1e5.mat'));
load(fullfile(dir, 'meas_1d_1e5.mat'));

% targets
ntrgs = 10000;
xtrgs = linspace(min(x), max(x), ntrgs);
opts.only_trgs = 1;

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
ker = SE_ker(dim,l);

% get accurate solution
opts.tol = 1e-12;
[y_true0, ytrg_true0, info0] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
fprintf('\n');

opts.tol = 1e-14;
[y_true, ytrg_true, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
fprintf('\n');
fprintf('max dd: %g\n', max(abs(ytrg_true0.mean - ytrg_true.mean)));
fprintf('rms dd: %g\n', rms(ytrg_true0.mean - ytrg_true.mean));
%save(fullfile(dir, 'se_1d_true.mat'), 'ytrg_true');
%save(fullfile(dir, 'se_1d_true0.mat'), 'ytrg_true0');

% EFGP
nns = 10;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.tol = 1e-2 * 10^(-i);
    opts.l2scaled = true;
    [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

efgp_1d.ts = ts;
efgp_1d.rms_errs = rms_errs;
efgp_1d.linf_errs = linf_errs;
save(fullfile(dir, 'efgp_1d.mat'), 'efgp_1d');




% SKI
fprintf('---SKI---\n');
nns = 4;
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

ski_1d.ts = ts;
ski_1d.rms_errs = rms_errs;
ski_1d.linf_errs = linf_errs;
save(fullfile(dir, 'ski_1d.mat'), 'ski_1d');





% FLAM
fprintf('---FLAM---\n');
nns = 5;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
opts.v = true;
for i=1:nns
    opts.tol = 1e-4 * 10^(-i);
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_1d.ts = ts;
flam_1d.rms_errs = rms_errs;
flam_1d.linf_errs = linf_errs;
save(fullfile(dir, 'flam_1d.mat'), 'flam_1d');







%RLCM
fprintf('---RLCM---\n');
nns = 4;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.rank = 30 * i;
    [y, ytrg, info] = RLCM(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

rlcm_1d.ts = ts;
rlcm_1d.rms_errs = rms_errs;
rlcm_1d.linf_errs = linf_errs;
save(fullfile(dir, 'rlcm_1d.mat'), 'rlcm_1d');







% plotting
load(fullfile(dir, 'efgp_1d.mat'));
load(fullfile(dir, 'ski_1d.mat'));
load(fullfile(dir, 'flam_1d.mat'));
load(fullfile(dir, 'rlcm_1d.mat'));
hold on;
plot(log10(efgp_1d.rms_errs), log10(efgp_1d.ts), '-o');
plot(log10(ski_1d.rms_errs), log10(ski_1d.ts), '-o');
plot(log10(flam_1d.rms_errs), log10(flam_1d.ts), '-o');
plot(log10(rlcm_1d.rms_errs), log10(rlcm_1d.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;


% print for paper
fprintf('EFGP\n')
print_tikz(log10(efgp_1d.rms_errs), log10(efgp_1d.ts))
fprintf('\nSKI\n')
print_tikz(log10(ski_1d.rms_errs), log10(ski_1d.ts))
fprintf('\nFLAM\n')
print_tikz(log10(flam_1d.rms_errs), log10(flam_1d.ts))
fprintf('\nRLCM\n')
print_tikz(log10(rlcm_1d.rms_errs), log10(rlcm_1d.ts))
