rng(1);
addpath ~/git/finufft/matlab/
run ../../../startup.m
% set directory for loading data and saving results
dir = "~/ceph/gp-shootout/big_example/data";
dir2 = '~/ceph/gp-shootout/big_example/results_oct14_2022/';

Nvals = [3e6; 1e7; 3e7; 1e8; 1e9];
sigvals = 0.1*sqrt(Nvals/3e6);

iNvals = 5;
isig =5;
iNtrg = 3;
alg = 'EFGP';

itol = 5;
tol = 10^(-itol);

fname = [dir2 alg '_iNvals' int2str(iNvals) '_isig' int2str(isig) '_iNtrg' int2str(iNtrg) '_itol' int2str(itol) '.mat'];

wavevec = [4; 3];
phase = 1.3;
f = @(x) cos(2*pi*x'* wavevec + phase);


% set tolerances and make sure we only evaluate posterior mean at trgs
opts2d.tol = tol;
opts2d.l2scaled = true;

% load data
N0 = Nvals(iNvals);
sigmatrue = sigvals(isig);

fstr = ['meas_2d_isig' int2str(isig) '_aug23_2022.bin'];
fid = fopen(fullfile(dir,fstr),'r');
tic,
meas = fread(fid,[N0,1],'double');
toc;
fclose(fid);

fid = fopen(fullfile(dir,'x_2d_aug23_2022.bin'),'r');
tic,
x = fread(fid,[2,N0],'double');
toc;
fclose(fid);




% subsample data
N = N0;
% meas = meas0(1:N);
% x = x0(:,1:N);

% gp regression
dim = 2;
l = 0.1;
nu = 1.5;
var = 1.0;
ker = Matern_ker(dim, nu, l, var);
%ker = SE_ker(dim, l, var);
ntrgs_per_d = 10^(iNtrg);
xtrgs = equispaced_grid(dim, ntrgs_per_d);

sigmasq = sigmatrue^2;

alpha =0;
memorygraph('start');
if(strcmpi(alg,'FLAMGP'))
    [y, ytrgs, info,alpha] = FLAMGP(x, meas, sigmasq, ker, xtrgs, opts2d);
elseif(strcmpi(alg,'EFGP'))
    [y, ytrgs, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts2d);
end
[bytes est_times cpu_times cpu_usages labelstrings labeltimes] = ...
   memorygraph('get');
memorygraph('done');
mginfo = [];
mginfo.bytes = bytes;
mginfo.est_times = est_times;
mginfo.cpu_time = cpu_times;
mginfo.labelstrings = labelstrings;
mginfo.labeltimes = labeltimes;
truemean = f(xtrgs);
err_sig = rms(truemean(:)-ytrgs.mean(:));

save(fname, 'l','dim', 'N','nu','var','sigmatrue','opts2d', ...
   'ytrgs','ntrgs_per_d','info', 'mginfo',...
   'alpha', 'err_sig','-v7.3');
%exit


