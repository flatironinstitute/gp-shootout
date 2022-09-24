rng(1);
sigmatrue = 0.3;  % used to regress

% set directory for loading data and saving results
dir = "~/gp-shootout/results/philip/efgp_tables/data";
    %%%%%%%%%%%%%%save(fullfile(dir, filename), 'err_rms')

save(fullfile(dir, 'sigmatrue.mat'),'sigmatrue');

% generate data (x, meas) for 1e7 points in 1, 2, and 3 dimensions
N = 1e7;

% 1d, N=1e7
dim = 1;
wavevec = 3;
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_1d_1e7.mat'), 'meas');
save(fullfile(dir, 'truemeas_1d_1e7.mat'), 'truemeas');
save(fullfile(dir, 'x_1d_1e7.mat'), 'x');
save(fullfile(dir, 'wavevec_1d_1e7.mat'),'wavevec');

% 2d, N = 1e7
dim = 2;
wavevec = [4; 3];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_2d_1e7.mat'), 'meas');
save(fullfile(dir, 'truemeas_2d_1e7.mat'), 'truemeas');
save(fullfile(dir, 'x_2d_1e7.mat'), 'x');
save(fullfile(dir, 'wavevec_2d_1e7.mat'),'wavevec');

% 3d, N = 1e7
dim = 3;
wavevec = [3; 7; 2];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_3d_1e7.mat'), 'meas');
save(fullfile(dir, 'truemeas_3d_1e7.mat'), 'truemeas');
save(fullfile(dir, 'x_3d_1e7.mat'), 'x');
save(fullfile(dir, 'wavevec_3d_1e7.mat'),'wavevec');





% 1d, N=1e8
[x, meas, truemeas] = get_randdata(dim, 1e8, f, sigmatrue);
save(fullfile(dir, 'meas_1d_1e8.mat'), 'meas');
save(fullfile(dir, 'x_1d_1e8.mat'), 'x');

% 2d, N = 1e8
[x, meas, truemeas] = get_randdata(dim, 1e8, f, sigmatrue);
save(fullfile(dir, 'meas_2d_1e8.mat'), 'meas');
save(fullfile(dir, 'x_2d_1e8.mat'), 'x');

% 3d, N = 1e8
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_3d_1e8.mat'), 'meas');
save(fullfile(dir, 'x_3d_1e8.mat'), 'x', '-v7.3'); % we need the -v7.3 otherwise we get this error message: Warning: Variable 'x' was not saved. For variables larger than 2GB use MAT-file version 7.3 or later. 


