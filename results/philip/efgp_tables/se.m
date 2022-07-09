% numerical results for efgp using squared-exponential kernel in 1, 2, 3
% dimensions. the results are used to construct tables in the paper. 
rng(1);

% sigma used to generate data and to be used for regression
dir = "~/gp-shootout/results/philip/efgp_tables/data";
load(fullfile(dir, 'sigmatrue.mat'));


% 1d
% load data
fprintf("\n dim=1 \n");
load(fullfile(dir, 'x_1d_1e7.mat'));
load(fullfile(dir, 'meas_1d_1e7.mat'));
dim = 1;
l = 0.1;
ker = SE_ker(dim,l);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);
opts.only_trgs = 1;

nns = 7;
for i=1:nns
    N = 10^i;
    % subsample
    x_i = x(1e7 * (1:N)/N);
    meas_i = meas(1e7 * (1:N)/N);
    
    sigmasq = sigmatrue^2;
    opts.tol = 1e-8;
    [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    
    % reference calculation
    opts.tol = 1e-10;
    [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);

    %disp(size(info.xis));
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
    m = (numel(info.xis) -  1) / 2;
    %%%disp(m);
%         fprintf('EFGP rms at targets %.3g, time: %.3g\n', rms(ytrgs.mean-ytrgs_true.mean), info.cpu_time.total);
%         fprintf('rms  at targets %.3g\n', rms(ytrgs.mean-ytrgs_true.mean));
%         fprintf('linf at targets %.3g\n', max(abs(ytrgs.mean-ytrgs_true.mean)));
%         fprintf('total time   %.3g\n', info.cpu_time.total);
%         fprintf('precomp time %.3g\n', info.cpu_time.precomp);
%         fprintf('cg time      %.3g\n', info.cpu_time.cg);
%         fprintf('mean time    %.3g\n', info.cpu_time.mean);
%         fprintf('mean/target  %.3g\n', info.cpu_time.mean / (ntrgs_per_d^dim));

     fprintf("$ 10^{%d}$ & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
         log10(N), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
     filename = sprintf('se_%gd_info_1e%g.mat', dim, log10(N));
     save(fullfile(dir, filename), 'info')
     filename = sprintf('se_%gd_rms_err_1e%g.mat', dim, log10(N));
     save(fullfile(dir, filename), 'err_rms')

end


% 2d
% load data
fprintf("\n dim=2 \n");
load(fullfile(dir, 'x_2d_1e7.mat'));
load(fullfile(dir, 'meas_2d_1e7.mat'));
dim = 2;
l = 0.1;
ker = SE_ker(dim,l);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);

nns = 7;
for i=1:nns
    N = 10^i;

    % subsample
    x_i = x(:, 1e7 * (1:N)/N);
    meas_i = meas(1e7 * (1:N)/N);
    x_i(:, 1) = x(:, 1);
    meas_i(1) = meas(1);

    % regression
    sigmasq = sigmatrue^2;
    clear opts;
    opts.tol = 1e-8;
    
    [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    opts.tol = 1e-10;
    [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
    m = (numel(info.xis) -  1) / 2;
    
    %fprintf("$ 10^{%d}$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
    %    log10(N), info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    fprintf("$ 10^{%d}$ & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
         log10(N), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    filename = sprintf('se_%gd_info_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'info');
    filename = sprintf('se_%gd_rms_err_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'err_rms');
    end


% 3d
% load data
fprintf("\n dim=3 \n");
dim = 3;
load(fullfile(dir, 'x_3d_1e7.mat'));
load(fullfile(dir, 'meas_3d_1e7.mat'));
l = 0.1;
ker = SE_ker(dim,l);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);

nns = 7;
for i=1:nns
    N = 10^i;
    % subsample
    x_i = x(:, 1e7 * (1:N)/N);
    meas_i = meas(1e7 * (1:N)/N);

    sigmasq = sigmatrue^2;
    opts.tol = 1e-4;
    [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    m = (numel(info.xis) -  1) / 2;
    
    % reference calculation
    opts.tol = 1e-6;
    [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
    
    fprintf("$ 10^{%d}$ & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ \\\\ \n", ...
         log10(N), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    filename = sprintf('se_%gd_info_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'info');
    filename = sprintf('se_%gd_rms_err_1e%g.mat', dim, log10(N));
    save(fullfile(dir, filename), 'err_rms');
end



% print stuff
for dim = 1:3
    fprintf('\ndim=%g\n',dim)
    for i=1:7
        N = 10^(i);
        % load info 
        filename = sprintf('se_%gd_info_1e%g.mat', dim, i);
        load(fullfile(dir, filename));
        % load errors
        filename = sprintf('se_%gd_rms_err_1e%g.mat', dim, log10(N));
        load(fullfile(dir, filename));
        % print for table
        fprintf("$ 10^{%d}$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %d $ \\\\ \n", ...
            log10(N), info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_rms);
    end
end
