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

% % subsample for testing
% N = 100;
% meas = meas(1:N);
% x = x(:, 1:N);

% targets
ntrgs = 30;
xtrgs = equispaced_grid(dim, ntrgs);
opts.only_trgs = 1;

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
var = 1;
nu = 0.5;
ker = Matern_ker(dim, nu, l, var);

% Accurate solution via efgp
opts.tol = 0.001;
[y_true0, ytrg_true0, info0] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
save(fullfile(dir, 'matern_3d_true.mat'), 'ytrg_true0');
fprintf('finished ground truth 1\n');

opts.tol = 0.0005;
[y_true, ytrg_true, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
fprintf('finished ground truth 2\n');
save(fullfile(dir, 'matern_3d_true.mat'), 'ytrg_true');

% % Accurate solution via flam
% opts.no_proxy = true;
% opts.tol = 1e-3; 
% opts.v = true;
% [y_true0, ytrg_true0, info0] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
% fprintf('finished ground truth 1\n');
% save(fullfile(dir, 'matern_3d_true0.mat'), 'ytrg_true0');
% 
% opts.tol = 1e-4;
% [y_true, ytrg_true, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
% fprintf('finished ground truth 2\n');
% save(fullfile(dir, 'matern_3d_true.mat'), 'ytrg_true');

save(fullfile(dir, 'matern_3d_x.mat'), 'x');
save(fullfile(dir, 'matern_3d_meas.mat'), 'meas');

load(fullfile(dir, 'matern_3d_true.mat'));
load(fullfile(dir, 'matern_3d_true0.mat'));
fprintf('max dd: %g\n', max(abs(ytrg_true.mean - ytrg_true0.mean)));
fprintf('rmse: %g\n', rms(ytrg_true.mean - ytrg_true0.mean));


% EFGP
nns = 3;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    fprintf('efgp %d\n', i);
    opts.tol = 10^(-i);
    [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
    
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

efgp_3d_matern.ts = ts;
efgp_3d_matern.rms_errs = rms_errs;
efgp_3d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'efgp_3d_matern.mat'), 'efgp_3d_matern');


% SKI
nns = 3;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    fprintf('ski %d\n', i);
    opts.grid_size = 5^i;
    [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrgs, opts);

    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

ski_3d_matern.ts = ts;
ski_3d_matern.rms_errs = rms_errs;
ski_3d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'ski_3d_matern.mat'), 'ski_3d_matern');



% FLAM
nns = 2;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
opts.v = true;
for i=1:nns
    fprintf('flam %d\n', i);
    opts.tol = 10^(-2 - i); 
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_3d_matern.ts = ts;
flam_3d_matern.rms_errs = rms_errs;
flam_3d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'flam_3d_matern.mat'), 'flam_3d_matern');





%RLCM
nns = 4;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    fprintf('rlcm %d\n', i);
    opts.rank = 150 * i;
    [y, ytrg, info] = RLCM(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

rlcm_3d_matern.ts = ts;
rlcm_3d_matern.rms_errs = rms_errs;
rlcm_3d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'rlcm_3d_matern.mat'), 'rlcm_3d_matern');





% plotting
load(fullfile(dir, "efgp_3d_matern.mat"));
load(fullfile(dir, "ski_3d_matern.mat"));
load(fullfile(dir, "flam_3d_matern.mat"));
load(fullfile(dir, "rlcm_3d_matern.mat"));
hold on;
plot(log10(efgp_3d_matern.rms_errs), log10(efgp_3d_matern.ts), '-o');
plot(log10(ski_3d_matern.rms_errs), log10(ski_3d_matern.ts), '-o');
plot(log10(flam_3d_matern.rms_errs), log10(flam_3d_matern.ts), '-o');
plot(log10(rlcm_3d_matern.rms_errs), log10(rlcm_3d_matern.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;




% print for paper
fprintf('EFGP\n')
print_tikz(log10(efgp_3d_matern.rms_errs), log10(efgp_3d_matern.ts))
fprintf('\nSKI\n')
print_tikz(log10(ski_3d_matern.rms_errs), log10(ski_3d_matern.ts))
fprintf('\nFLAM\n')
print_tikz(log10(flam_3d_matern.rms_errs), log10(flam_3d_matern.ts))
fprintf('\nRLCM\n')
print_tikz(log10(rlcm_3d_matern.rms_errs), log10(rlcm_3d_matern.ts))