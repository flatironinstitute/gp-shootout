rng(1);
sigmatrue = 0.3;

% set directory for saving data
dir = "~/ceph/gp-shootout/big_example/data";
save(fullfile(dir, 'sigmatrue.mat'), 'sigmatrue');

% 2d
N = 1e9;
dim = 2;
wavevec = [4; 3];
f = @(x) cos(2*pi*x'* wavevec + 1.3);
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
% save(fullfile(dir, 'meas_2d.mat'), 'meas', '-v7.3');
% save(fullfile(dir, 'x_2d.mat'), 'x', '-v7.3');

tic,
fid = fopen(fullfile(dir,'meas_2d.bin'),'w');
fwrite(fid,meas,'double');
toc;
fclose(fid)

tic,
fid = fopen(fullfile(dir,'x_2d.bin'),'w');
fwrite(fid,x,'double');
toc;
fclose(fid);

