rng(1);

% set directory for loading data and saving results
dir = [fileparts(mfilename('fullpath')) '/data'];
system(['mkdir -p ' dir]);
sigmatrue = 0.3;  % used to regress
save(fullfile(dir, 'sigmatrue.mat'),'sigmatrue');

% generate data (x, meas) for 1e7 points in 1, 2, and 3 dimensions
N = 1e7;

% 1d, N=1e7
dim = 1;
wavevec = 3.0;
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_1d_1e7.mat'), 'meas');
save(fullfile(dir, 'truemeas_1d_1e7.mat'), 'truemeas');
save(fullfile(dir, 'x_1d_1e7.mat'), 'x');
save(fullfile(dir, 'wavevec_1d_1e7.mat'),'wavevec');
save(fullfile(dir, 'f_1d.mat'), 'f')

% 2d, N = 1e7
dim = 2;
wavevec = [1; 2];
wavevec = 3.0 * wavevec / norm(wavevec);
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_2d_1e7.mat'), 'meas');
save(fullfile(dir, 'truemeas_2d_1e7.mat'), 'truemeas');
save(fullfile(dir, 'x_2d_1e7.mat'), 'x');
save(fullfile(dir, 'wavevec_2d_1e7.mat'),'wavevec');
save(fullfile(dir, 'f_2d.mat'), 'f')

% 3d, N = 1e7
dim = 3;
wavevec = [1; 3; 2];
wavevec = 3.0 * wavevec / norm(wavevec);
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_3d_1e7.mat'), 'meas');
save(fullfile(dir, 'truemeas_3d_1e7.mat'), 'truemeas');
save(fullfile(dir, 'x_3d_1e7.mat'), 'x');
save(fullfile(dir, 'wavevec_3d_1e7.mat'),'wavevec');
save(fullfile(dir, 'f_3d.mat'), 'f')