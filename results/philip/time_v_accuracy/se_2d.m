% accuracy vs time for for fixed N

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
ker = SE_ker(dim,l);

% get accurate solution
opts.tol = 1e-14;
[y, ytrg_true, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
opts.tol = 1e-15;
[y, ytrg_true2, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
%%%scatter(x, meas); hold on; plot(xtrgs, ytrg_true.mean); hold off
%%%[y, ytrg_true, ~] = naive_gp(x, meas, sigmasq, ker, xtrgs, []);
fprintf('max dd: %g\n', max(abs(ytrg_true.mean - ytrg_true2.mean)));



% EFGP
nns = 10;
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
save('efgp_2d.mat','efgp_2d');


% SKI
nns = 4;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    opts.grid_size = 5^i;
    [y, ytrg, info] = SKI(x, meas, sigmasq, ker, xtrgs, opts);

    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

ski_2d.ts = ts;
ski_2d.rms_errs = rms_errs;
ski_2d.linf_errs = linf_errs;
save('ski_2d.mat','ski_2d');



% FLAM
nns = 1;
ts = zeros(nns, 1);
linf_errs = zeros(nns, 1);
rms_errs = zeros(nns, 1);
for i=1:nns
    [y, ytrg, info] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts);
    ts(i) = info.cpu_time.total;
    rms_errs(i) = rms(ytrg.mean - ytrg_true.mean);
    linf_errs(i) = max(abs(ytrg.mean - ytrg_true.mean));
end

flam_2d.ts = ts;
flam_2d.rms_errs = rms_errs;
flam_2d.linf_errs = linf_errs;
save('flam_2d.mat','flam_2d');



% plotting
load("efgp_2d.mat");
load("ski_2d.mat");
load("flam_2d.mat");
hold on;
plot(log10(efgp_2d.rms_errs), log10(efgp_2d.ts), '-o');
plot(log10(ski_2d.rms_errs), log10(ski_2d.ts), '-o');
plot(log10(flam_2d.rms_errs), log10(flam_2d.ts), '-o');
set ( gca, 'xdir', 'reverse' );
hold off;
