% accuracy vs time for for fixed N
clear opts;

% set directory for saving results and loading data
dir = [fileparts(mfilename('fullpath')) '/data'];

% sigma used to generate data and to be used for regression
load(fullfile(dir, 'sigmatrue.mat'));

% load data
load(fullfile(dir, 'x_1d_1e5.mat'));
load(fullfile(dir, 'meas_1d_1e5.mat'));

% targets
ntrgs = 10000;
xtrgs = linspace(0, 1, ntrgs);
opts.only_trgs = 1;

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
dim = 1;
var = 1;
nu = 0.5;
ker = Matern_ker(dim, nu, l, var);

% % get accurate solution
% opts.tol = 1e-14;
% [y0, ytrg_true0, info0] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
% opts.tol = 1e-15;
% [y, ytrg_true, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
% save(fullfile(dir, 'matern_1d_true.mat'), 'ytrg_true');
% save(fullfile(dir, 'matern_1d_true0.mat'), 'ytrg_true0');
% save(fullfile(dir, 'matern_1d_x.mat'), 'x');
% save(fullfile(dir, 'matern_1d_meas.mat'), 'meas');

% load reference solution and data
load(fullfile(dir, 'matern_1d_true.mat'));
load(fullfile(dir, 'matern_1d_true0.mat'));
load(fullfile(dir, 'matern_1d_x.mat'));
load(fullfile(dir, 'matern_1d_meas.mat'));
fprintf('max dd: %g\n', max(abs(ytrg_true.mean - ytrg_true0.mean)));


% EFGP
nns = 5;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.l2scaled = true;
    opts.tol = 1e-2 * 10^(-i);
    [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
    
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

efgp_1d_matern.ts = ts;
efgp_1d_matern.rms_errs = rms_errs;
efgp_1d_matern.linf_errs = linf_errs;
save('efgp_1d_matern.mat','efgp_1d_matern');

return


% SKI
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

ski_1d_matern.ts = ts;
ski_1d_matern.rms_errs = rms_errs;
ski_1d_matern.linf_errs = linf_errs;
save('ski_1d_matern.mat','ski_1d_matern');



% FLAM
nns = 6;
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

flam_1d_matern.ts = ts;
flam_1d_matern.rms_errs = rms_errs;
flam_1d_matern.linf_errs = linf_errs;
save('flam_1d_matern.mat','flam_1d_matern');

% RLCM
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
rlcm_1d_matern.ts = ts;
rlcm_1d_matern.rms_errs = rms_errs;
rlcm_1d_matern.linf_errs = linf_errs;
save(fullfile(dir, 'rlcm_1d_matern.mat'), 'rlcm_1d_matern');


% plotting
load("efgp_1d_matern.mat");
load("ski_1d_matern.mat");
load("flam_1d_matern.mat");
load("rlcm_1d_matern.mat");
hold on;
plot(log10(efgp_1d_matern.rms_errs), log10(efgp_1d_matern.ts), '-o');
plot(log10(ski_1d_matern.rms_errs), log10(ski_1d_matern.ts), '-o');
plot(log10(flam_1d_matern.rms_errs), log10(flam_1d_matern.ts), '-o');
plot(log10(rlcm_1d_matern.rms_errs), log10(rlcm_1d_matern.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;




% print for paper
fprintf('EFGP\n')
print_tikz(log10(efgp_1d_matern.rms_errs), log10(efgp_1d_matern.ts))
fprintf('\nSKI\n')
print_tikz(log10(ski_1d_matern.rms_errs), log10(ski_1d_matern.ts))
fprintf('\nFLAM\n')
print_tikz(log10(flam_1d_matern.rms_errs), log10(flam_1d_matern.ts))
fprintf('\nRLCM\n')
print_tikz(log10(rlcm_1d_matern.rms_errs), log10(rlcm_1d_matern.ts))



