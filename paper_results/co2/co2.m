% pull in data
[x0, meas0] = get_CO2data(2000000);

lat_max = 65;
lat_min = -50;
mask = (x0(2, :) < lat_max) & (x0(2, :) > lat_min);
x0 = x0(:, mask);
meas0 = meas0(mask);

% demean measurements
meas_cen = mean(meas0);
meas = meas0 - meas_cen;

dim = 2;
xtrgs_x = linspace(-180, 180, 300);
xtrgs_y = linspace(lat_min, lat_max, 300);
[xx, yy] = meshgrid(xtrgs_x, xtrgs_y);
xtrgs = [reshape(xx, 1, []); reshape(yy, 1, [])];

l = 50;
sigmasq = 1.0;
ker = SE_ker(dim, l);
ker.k = @(x) 3^2 * ker.k(x);
ker.khat = @(x) 3^2 * ker.khat(x);
    
clear opts
opts.tol = 1e-7;
opts.l2scaled = true;
[y1, ytrg1, info1] = EFGP(x0, meas, sigmasq, ker, xtrgs, opts);
fprintf('finished run 1\n')

opts.tol = 1e-8;
[y0, ytrg0, info0] = EFGP(x0, meas, sigmasq, ker, xtrgs, opts);
fprintf('finished run true\n')
   
rms_err1 = rms(ytrg1.mean - ytrg0.mean) / rms(ytrg0.mean);

% add back mean
ytrg1.mean = ytrg1.mean + meas_cen;
ytrg0.mean = ytrg0.mean + meas_cen;

% plot
geolimits([-70 70], [-180 180])

subplot(1,3,2)
geoscatter(xtrgs(2,:), xtrgs(1,:), 10, ytrg1.mean, 'filled');
caxis manual
caxis([min(meas0, [], 'all') max(meas0, [], 'all')]);
colorbar;

subplot(1,3,1)
geoscatter(x0(2,:), x0(1,:), 1, meas0, '^');
caxis manual
caxis([min(meas0, [], 'all') max(meas0, [], 'all')]);
colorbar;

subplot(1,3,3)
geoscatter(xtrgs(2,:), xtrgs(1,:), 10, ytrg1.mean - ytrg0.mean, '^');
colorbar;


hold off;


fprintf('N: %g\n', numel(meas));

fprintf('l: %g\n', l);
fprintf('total time1: %g\n', info1.cpu_time.total);
fprintf('rms error1: %g\n', rms_err1);
fprintf('iterations: %d\n', info1.iter);