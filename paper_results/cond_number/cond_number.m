rng(1);

% sigma used to generate data and to be used for regression
%dir = "~/gp-shootout/paper_results/efgp_tables/data";
dir = [fileparts(mfilename('fullpath')) '/../efgp_tables/data'];
load(fullfile(dir, 'sigmatrue.mat'));

% load data
load(fullfile(dir, 'x_1d_1e7.mat'));
load(fullfile(dir, 'meas_1d_1e7.mat'));

% parameters
dim = 1;
l = 0.1;
ker = SE_ker(dim,l);
ntrgs_per_d = 100;
xtrgs = equispaced_grid(dim, ntrgs_per_d);
opts.only_trgs = 1;
opts.dense = true;

% condition number as a function of number of data points
nns = 6;
kappas = zeros(nns, 1);
bounds = zeros(nns, 1);
trues = zeros(nns, 1);
for i=1:nns
    N = 10^i;
    fprintf('N: %d\n', N)
    % subsample
    x_i = x(1e7 * (1:N)/N);
    meas_i = meas(1e7 * (1:N)/N);
    
    sigmasq = sigmatrue^2;
    opts.tol = 1e-10;
    [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    
    % reference calculation
    opts.tol = 1e-10;
    [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);

    % accuracy
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));

    m = numel(info.xis);
    kappas(i) = log10(cond(info.A + sigmasq * eye(m)));
    bounds(i) = log10(N / sigmasq + 1 + sigmasq/2);

    % check exact covariance matrix for small N
    if N <= 1e3 % occasional performance issues for 1e4
        K = densekermat(ker.k,x_i);
        trues(i) = log10(cond(K + sigmasq * eye(N)));
    end
end


% print for tikzpicture
fprintf('k(X^*X + sigma^2*I):\n')
print_tikz(1:1:numel(kappas), kappas)
fprintf('upper bounds:\n')
print_tikz(1:1:numel(bounds), bounds)
fprintf('k(K + sigma^2I):\n')
print_tikz(1:1:numel(trues), trues)






% heatmaps with condition number as function of N and sigma^2

% parameters
nns = 4;
nsigs = 11;
kappas = zeros(nns, nsigs);
bounds = zeros(nns, nsigs);
ratios = zeros(nns, nsigs);
ns = zeros(nns, 1);
sigs = linspace(0.1, 2.1, nsigs);
for i=1:nns
    N = 10^(i+1);
    ns(i) = log10(N);
    % subsample
    x_i = x(1e7 * (1:N)/N);
    meas_i = meas(1e7 * (1:N)/N);
    
    for j=1:nsigs
        sigmasq = sigs(j)^2;
        fprintf('N: %d, sigma^2: %0.2f\n', N, sigmasq);
        opts.tol = 1e-8;
        [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
        
        % reference calculation
        opts.tol = 1e-10;
        [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    
        % check accuracy
        err_rms = rms(ytrgs.mean-ytrgs_true.mean);
        err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
    
        m = numel(info.xis);
        kappas(i, j) = log10(cond(info.A + sigmasq * eye(m)));
        bounds(i, j) = log10(N / sigmasq + 1 + sigmasq/2);
        ratios(i, j) = (exp(kappas(i, j)) ./ exp(bounds(i, j)));
    end    
end


% plot heatmaps
subplot(1, 3, 1);
imagesc(kappas);
set(gca, 'XTick', 1:nsigs, 'XTickLabel', sigs)
set(gca, 'YTick', 1:nns,   'YTickLabel', ns)
set(gca,'YDir','normal')
caxis manual
caxis([0 max(kappas, [], 'all')]);
colorbar;

subplot(1, 3, 2);
imagesc(bounds);
set(gca, 'XTick', 1:nsigs, 'XTickLabel', sigs)
set(gca, 'YTick', 1:nns,   'YTickLabel', ns)
set(gca,'YDir','normal')
caxis manual
caxis([0 max(kappas, [], 'all')]);
colorbar;

subplot(1, 3, 3);
imagesc(ratios);
set(gca, 'XTick', 1:nsigs, 'XTickLabel', sigs)
set(gca, 'YTick', 1:nns,   'YTickLabel', ns)
set(gca,'YDir','normal')
colorMap = [linspace(0.2, 0.7,256)', linspace(0.2, 0.7,256)', ones(256, 1)];
colormap(gca, colorMap);
colorbar;