rng(1);

% sigma used to generate data and to be used for regression
dir = "~/gp-shootout/results/philip/efgp_tables/data";
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
nns = 2;
kappas = zeros(nns, 1);
bounds = zeros(nns, 1);
trues = zeros(nns, 1);
for i=1:nns
    fprintf('%d\n', i)
    N = 10^i;
    % subsample
    x_i = x(1e7 * (1:N)/N);
    meas_i = meas(1e7 * (1:N)/N);
    
    sigmasq = sigmatrue^2;
    opts.tol = 1e-10;
    [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    
    % reference calculation
    opts.tol = 1e-10;
    [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);

    %disp(size(info.xis));
    err_rms = rms(ytrgs.mean-ytrgs_true.mean);
    err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));

    m = numel(info.xis);
    X = info.X;
    XtX = X' * X;
    kappas(i) = log10(cond(info.A + sigmasq * eye(m)));
    bounds(i) = log10(N / sigmasq + 1 + sigmasq/2);

    % check exact covariance matrix 
    K = densekermat(ker.k,x_i);
    trues(i) = log10(cond(K + sigmasq * eye(N)));
    
end

return

% print for tikzpicture
fprintf('kappas:\n')
print_tikz(1:1:numel(kappas), kappas)
fprintf('bounds:\n')
print_tikz(1:1:numel(bounds), bounds)
fprintf('trues:\n')
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
        fprintf('%d, %d\n', i, j);
        sigmasq = sigs(j)^2;
        opts.tol = 1e-8;
        [y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
        
        % reference calculation
        opts.tol = 1e-10;
        [y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);
    
        %disp(size(info.xis));
        err_rms = rms(ytrgs.mean-ytrgs_true.mean);
        err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));
    
        m = numel(info.xis);
        X = info.X;
        XtX = X' * X;
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
%%xtickformat('%.2f')
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
%%%% log of ratio, not ratio of logs
%%%imagesc(kappas ./ bounds);
imagesc(ratios);
set(gca, 'XTick', 1:nsigs, 'XTickLabel', sigs)
set(gca, 'YTick', 1:nns,   'YTickLabel', ns)
set(gca,'YDir','normal')
%caxis manual
%caxis(log([min(ratios, [] , 'all') max(ratios, [], 'all')]));
%colormap(gca, winter);
colorMap = [linspace(0.2, 0.7,256)', linspace(0.2, 0.7,256)', ones(256, 1)];
colormap(gca, colorMap);
colorbar;




% parameters
N = 10^2;
x_i = x(1e7 * (1:N)/N);
meas_i = meas(1e7 * (1:N)/N);
sig = 2.0;
sigmasq = sig^2;
opts.tol = 1e-8;
[y, ytrgs, info] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);

% reference calculation
opts.tol = 1e-10;
[y_true, ytrgs_true, info_true] = EFGP(x_i, meas_i, sigmasq, ker, xtrgs, opts);

%disp(size(info.xis));
err_rms = rms(ytrgs.mean-ytrgs_true.mean);
err_linf = max(abs(ytrgs.mean-ytrgs_true.mean));

m = numel(info.xis);
X = info.X;
XtX = X' * X;
kappa = log10(cond(info.A + sigmasq * eye(m)));
bound = log10(N / sigmasq + 1 + sigmasq/2);
kappa
bound

