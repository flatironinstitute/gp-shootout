% numerical results for efgp using squared-exponential kernel in 1, 2, 3
% dimensions. the results are used to construct tables in the paper. 
clear
rng(1);

% sigma used to generate data and to be used for regression
dir = [fileparts(mfilename('fullpath')) '/data'];
load(fullfile(dir, 'sigmatrue.mat'));

% set opts
clear opts
clear opts_ref
opts.l2scaled = 1;
opts_ref.l2scaled = 1;
opts_ref.only_trgs = 1;

for dim=1:3
    fprintf("\n dim=%d \n", dim);
    filename = sprintf('x_%gd_1e7.mat', dim);
    load(fullfile(dir, filename));
    filename = sprintf('meas_%gd_1e7.mat', dim);
    load(fullfile(dir, filename));
    filename = sprintf('truemeas_%gd_1e7.mat', dim);
    load(fullfile(dir, filename));
    filename = sprintf('wavevec_%gd_1e7.mat', dim);
    load(fullfile(dir, filename),'wavevec');
    filename = sprintf('f_%gd.mat', dim);
    load(fullfile(dir, filename));

    nu = 0.5;
    l = 0.1;  % alex added since missing
    var = 1.0;
    ker = Matern_ker(dim, nu, l, var);
    ntrgs_per_d = 60;
    xtrgs = equispaced_grid(dim, ntrgs_per_d);
    ntrgs = ntrgs_per_d^(dim);
    ftrgs = f(xtrgs);
    % generate noisy samples at targets
    ytrgs_test = f(xtrgs) + sigmatrue * randn(ntrgs, 1);
    
    if dim == 1
        opts.tol = 1e-4;
    elseif dim == 2
        opts.tol = 1e-3;
    else
        opts.tol = 0.5*1e-2;   % since 1e-3 was taking >> 2 hrs for 3d :(
    end
    opts_ref.tol = opts.tol/10;

    for i=3:7
        N = 10^i;
        % subsample
        inds = numel(meas) * (1:N)/N;
        x_i = x(:,inds);
        meas_i = meas(inds);
        truemeas_i = truemeas(inds);
        sigmasq = sigmatrue^2;
        [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
        
        % reference calculation
        if N<=1e4
            [y_true, ytrgs_true, info_true] = naive_gp(x_i, meas_i, sigmasq, ker, xtrgs, opts_ref);
        else
            [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts_ref);
        end

        err_rmse_ex = rms(ytrgs_test - ytrgs_true.mean);
        err_eepm = rms(ytrgs.mean-ytrgs_true.mean);
        err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
        err_rms = rms(ytrgs_test - ytrgs.mean);
        m = (numel(info.xis) -  1) / 2;
        
        % print for latex table
        fprintf("$ 10^{%d}$ & $ %d $ & %s & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ & $ %.1d $ & $ %.1d $ \\\\ \n", ...
                log10(N), dim, 'Mat $1/2$', m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_eepm, err_rms, err_rmse_ex);
        info.opts = opts; info.opts_ref = opts_ref;    % crucial to save
        filename = sprintf('mat_%gd_info_1e%g.mat', dim, log10(N));
        save(fullfile(dir, filename), 'info')
        filename = sprintf('mat_%gd_eepm_err_1e%g.mat', dim, log10(N));
        save(fullfile(dir, filename), 'err_eepm');
        filename = sprintf('mat_%gd_rms_err_1e%g.mat', dim, log10(N));
        save(fullfile(dir, filename), 'err_rms');
        filename = sprintf('mat_%gd_err_rmse_ex_1e%g.mat', dim, log10(N));
        save(fullfile(dir, filename), 'err_rmse_ex');
    end

end



if 0   % print stuff ..... obsolete, see nicelatextable.m
for dim = 1:3
    for i=3:7
        N = 10^(i);
        % load info 
        filename = sprintf('mat_%gd_info_1e%g.mat', dim, log10(N));
        load(fullfile(dir, filename));
        filename = sprintf('mat_%gd_rms_err_1e%g.mat', dim, log10(N));
        load(fullfile(dir, filename));
        filename = sprintf('mat_%gd_eepm_err_1e%g.mat', dim, log10(N));
        load(fullfile(dir, filename));
        filename = sprintf('mat_%gd_err_rmse_ex_1e%g.mat', dim, log10(N));
        load(fullfile(dir, filename));

        m = (numel(info.xis) -  1) / 2;
        % print for table
        fprintf("$ 10^{%d}$ & $ %d $ & %s & $ %d $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %.1d $ & $ %.1d $ & $ %.1d $ \\\\ \n", ...
             log10(N), dim, 'Mat $1/2$', m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, err_eepm, err_rms, err_rmse_ex);
    end
    fprintf("\\hline \\hline \n")
end
end
