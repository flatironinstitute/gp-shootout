rng(1);
sigmatrue = 0.5;  % used to regress
save('sigmatrue.mat','sigmatrue');

% generate data (x, meas) for 1e7 points in 1, 2, and 3 dimensions
N = 1e5;

% 1d
dim = 1;
wavevec = 3;
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save('meas_1d_1e5.mat','meas');
save('x_1d_1e5.mat','x');

% 2d
dim = 2;
wavevec = [4; 3];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save('meas_2d_1e5.mat','meas');
save('x_2d_1e5.mat','x');
