rng(1);

% set directory for loading data and saving results
dir = "~/gp-shootout/results/philip/efgp_tables/data";
load(fullfile(dir, 'sigmatrue.mat'));


% set tolerances and make sure we only evaluate posterior mean at trgs
opts.tol = 1e-6;
opts.only_trgs = true;
opts_true.tol = 1e-7;
opts_true.only_trgs = true;

opts2d.tol = 1e-3;
opts2d.only_trgs = true;
opts2d_true.tol = 1e-4;
opts2d_true.only_trgs = true;

opts3d.tol = 1e-2;
opts3d.only_trgs = true;
opts3d_true.tol = 1e-3;
opts3d_true.only_trgs = true;



% 1d
% load data
load(fullfile(dir, 'x_1d_1e7.mat'));
load(fullfile(dir, 'meas_1d_1e7.mat'));
fprintf('-----1d-----\n');
dim = 1;
l = 0.1;
nu = 0.5;
var = 1.0;
ker = Matern_ker(dim, nu, l, var);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);

nns = 7;
for i=1:nns
    N = 10^i;
    % subsample
    xi = x(1e7 * (1:N)/N);
    measi = meas(1e7 * (1:N)/N);

    sigmasq = sigmatrue^2;
    [y, ytrgs, info] = EFGP(xi, measi, sigmasq, ker, xtrgs, opts);
    [ytrue, ytrgs_true, info_true] = EFGP(xi, measi, sigmasq, ker, xtrgs, opts_true);
    
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
%         fprintf('EFGP rms at targets %.3g, time: %.3g\n', rms(ytrgs.mean-ytrgs_true.mean), info.cpu_time.total);
%         fprintf('rms  at targets %.3g\n', rms(ytrgs.mean-ytrgs_true.mean));
%         fprintf('linf at targets %.3g\n', max(abs(ytrgs.mean-ytrgs_true.mean)));
%         fprintf('total time   %.3g\n', info.cpu_time.total);
%         fprintf('precomp time %.3g\n', info.cpu_time.precomp);
%         fprintf('cg time      %.3g\n', info.cpu_time.cg);
%         fprintf('mean time    %.3g\n', info.cpu_time.mean);
%         fprintf('mean/target  %.3g\n', info.cpu_time.mean / (ntrgs_per_d^dim));

    %fprintf("$ 10^{%d}$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %d $ \\\\ \n", ...
    %    log10(N), info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    m = (numel(info.xis) -  1) / 2;
    fprintf("$ 10^{%d}$ & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
         log10(N), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    

    filename = sprintf('mat_%gd_info_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'info')
    filename = sprintf('mat_%gd_rms_err_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'err_rms')
end





% 2d
% load data
load(fullfile(dir, 'x_2d_1e7.mat'));
load(fullfile(dir, 'meas_2d_1e7.mat'));
fprintf('-----2d-----\n');
dim = 2;
l = 0.1;
nu = 0.5;
var = 1.0;
ker = Matern_ker(dim, nu, l, var);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);
nns = 7;
for i=1:nns
    N = 10^i;
    % subsample
    xi = x(:, 1e7 * (1:N)/N);
    measi = meas(1e7 * (1:N)/N);

    sigmasq = sigmatrue^2;
    [y, ytrgs, info] = EFGP(xi, measi, sigmasq, ker, xtrgs, opts2d);
    [y_true, ytrgs_true, info_true] = EFGP(xi, measi, sigmasq, ker, xtrgs, opts2d_true);
    
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
    
    %fprintf("$ 10^{%d}$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %d $ \\\\ \n", ...
    %    log10(N), info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    m = (numel(info.xis) -  1) / 2;
    fprintf("$ 10^{%d}$ & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
         log10(N), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    

    filename = sprintf('mat_%gd_info_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'info')
    filename = sprintf('mat_%gd_rms_err_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'err_rms')
end








% 3d
% load data
load(fullfile(dir, 'x_3d_1e7.mat'));
load(fullfile(dir, 'meas_3d_1e7.mat'));
fprintf('-----3d-----\n');
dim = 3;
l = 0.1;
nu = 0.5;
var = 1.0;
ker = Matern_ker(dim, nu, l, var);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);

nns = 7;
for i=1:nns
    N = 10^i;
    % subsample
    xi = x(:, 1e7 * (1:N)/N);
    measi = meas(1e7 * (1:N)/N);

    sigmasq = sigmatrue^2;
    [y, ytrgs, info] = EFGP(xi, measi, sigmasq, ker, xtrgs, opts3d);
    opts3d_true.tol = opts3d.tol / 2.0;
    [ytrue, ytrgs_true, info_true] = EFGP(xi, measi, sigmasq, ker, xtrgs, opts3d_true);
    
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
    
    %%%fprintf("$ 10^{%d}$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %d $ \\\\ \n", ...
    %%%    log10(N), info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    m = (numel(info.xis) -  1) / 2;
    fprintf("$ 10^{%d}$ & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
         log10(N), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    
    
    filename = sprintf('mat_%gd_info_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'info')
    filename = sprintf('mat_%gd_rms_err_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'err_rms')
end




% print stuff
for dim = 1:3
    fprintf('\ndim=%g\n',dim)
    for i=1:7
        N = 10^(i);
        % load info 
        filename = sprintf('mat_%gd_info_1e%g.mat', dim, i);
        load(fullfile(dir, filename));
        % load errors
        filename = sprintf('mat_%gd_rms_err_1e%g.mat', dim, log10(N));
        load(fullfile(dir, filename));
        % print for table
        fprintf("$ 10^{%d}$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %d $ \\\\ \n", ...
            log10(N), info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    end
end
