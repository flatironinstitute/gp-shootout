% accuracy vs time for for fixed N
rng(1);

% data
N = 1e5;
f = @(x) cos(6*2*pi*x) / 2;
sigmatrue = 0.5;
dim = 2;
unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
wavevec = freqdata*unitvec;    % col vec
f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);% for ski
x(:, 1) = [0; 0];
x(:, N) = [1; 1];

% targets
ntrgs = 100;
xtrgs = equispaced_grid(dim, ntrgs);

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
var = 1;
nu = 0.5;
ker = Matern_ker(dim, nu, l, var);


% % get accurate solution
% opts.tol = 1e-4;
% [y0, ytrg_true0, info0] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
% opts.tol = opts.tol / 10;
% [y, ytrg_true, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
% fprintf('max dd: %g\n', max(abs(ytrg_true.mean - ytrg_true0.mean)));
% save('matern_2d_true.mat','ytrg_true');
% save('matern_2d_true0.mat','ytrg_true0');
% save('matern_2d_x.mat','x');
% save('matern_2d_meas.mat','meas');

load('matern_2d_true.mat');
load('matern_2d_x.mat');
load('matern_2d_meas.mat');


% EFGP
nns = 3;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.tol = 1e-1 * 10^(-i);
    [y, ytrg, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
    
    ts(i) = info.cpu_time(end);
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

efgp_2d_matern.ts = ts;
efgp_2d_matern.rms_errs = rms_errs;
efgp_2d_matern.linf_errs = linf_errs;
save('efgp_2d_matern.mat','efgp_2d_matern');


% SKI
nns = 2;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.grid_size = 5^i;
    [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrgs, opts);

    ts(i) = info.cpu_time(end);
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

ski_2d_matern.ts = ts;
ski_2d_matern.rms_errs = rms_errs;
ski_2d_matern.linf_errs = linf_errs;
save('ski_2d_matern.mat','ski_2d_matern');



% FLAM
nns = 1;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time(end);
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_2d_matern.ts = ts;
flam_2d_matern.rms_errs = rms_errs;
flam_2d_matern.linf_errs = linf_errs;
save('flam_2d_matern.mat','flam_2d_matern');



% plotting
load("efgp_2d_matern.mat");
load("ski_2d_matern.mat");
load("flam_2d_matern.mat");
hold on;
plot(log10(efgp_2d_matern.rms_errs), log10(efgp_2d_matern.ts), '-o');
plot(log10(ski_2d_matern.rms_errs), log10(ski_2d_matern.ts), '-o');
plot(log10(flam_2d_matern.rms_errs), log10(flam_2d_matern.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;
