rng(1);
sigmatrue = 0.3;  % used to regress
save('sigmatrue.mat','sigmatrue');

% generate data (x, meas) for 1e7 points in 1, 2, and 3 dimensions
N = 1e7;

% 1d
dim = 1;
wavevec = 3;
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save('meas_1d_1e7.mat','meas');
save('x_1d_1e7.mat','x');
save('wavevec_1d_1e7.mat','wavevec');


% 2d
dim = 2;
wavevec = [4; 3];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save('meas_2d_1e7.mat','meas');
save('x_2d_1e7.mat','x');
save('wavevec_2d_1e7.mat','wavevec');


% 3d
dim = 3;
wavevec = [3; 7; 2];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
save('meas_3d_1e7.mat','meas');
save('x_3d_1e7.mat','x');
save('wavevec_3d_1e7.mat','wavevec');

