rng(1);

% set directory for saving data
dir = "~/gp-shootout/results/philip/time_v_accuracy/data";

% sigma 
sigmatrue = 0.5;
save(fullfile(dir, 'sigmatrue.mat'), 'sigmatrue');


% generate data (x, meas) for 1e7 points in 1, 2, and 3 dimensions
N = 1e5;

% 1d
dim = 1;
wavevec = 3;
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_1d_1e5.mat'), 'meas');
save(fullfile(dir, 'x_1d_1e5.mat'), 'x');

% 2d
dim = 2;
wavevec = [4; 3];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_2d_1e5.mat'), 'meas');
save(fullfile(dir, 'x_2d_1e5.mat'), 'x');

% 3d
dim = 3;
wavevec = [2; 3; 5];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_3d_1e5.mat'), 'meas');
save(fullfile(dir, 'x_3d_1e5.mat'), 'x');
