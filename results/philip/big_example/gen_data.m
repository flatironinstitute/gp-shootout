rng(1);
sigmatrue = 0.3;

% set directory for saving data
dir = "~/gp-shootout/results/philip/big_example/data";
save(fullfile(dir, 'sigmatrue.mat'), 'sigmatrue');

% 2d
N = 5e8;
dim = 2;
wavevec = [4; 3];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save(fullfile(dir, 'meas_2d.mat'), 'meas', '-v7.3');
save(fullfile(dir, 'x_2d.mat'), 'x', '-v7.3');
