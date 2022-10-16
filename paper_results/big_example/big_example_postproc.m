clear
% addpath ~/git/finufft/matlab/
% run ../../../startup.m
% set directory for loading data and saving results
dir = '~/ceph/gp-shootout/big_example/results_sep5_2022/'; %directory where results are stored
dir2 = './results/'; % director for saving results
dir0 = "~/ceph/gp-shootout/big_example/data"; % directory where data is stored


Nvals = [3e6; 1e7; 3e7; 1e8; 1e9];
sigvals = 0.1*sqrt(Nvals/3e6);

wavevec = [4; 3];
phase = 1.3;
f = @(x) cos(2*pi*x'* wavevec + phase);


iNvals = 5;
isig = 1;
iNtrg = 3;

% load data
N0 = Nvals(iNvals);
sigmatrue = sigvals(isig);

sigmasq = sigmatrue^2;


fstr = ['meas_2d_isig' int2str(isig) '_aug23_2022.bin'];
fid = fopen(fullfile(dir0,fstr),'r');
tic,
meas = fread(fid,[N0,1],'double');
toc;
fclose(fid);

alg = 'FLAMGP';
alg2 = 'EFGP';
itol = 5;
fname_tol5 = [dir alg '_iNvals' int2str(iNvals) '_isig' int2str(isig) '_iNtrg' int2str(iNtrg) '_itol' int2str(itol) '.mat'];

itol = 7;
fname_tol7 = [dir alg '_iNvals' int2str(iNvals) '_isig' int2str(isig) '_iNtrg' int2str(iNtrg) '_itol' int2str(itol) '.mat'];

itol = 5;
fname2_tol5 = [dir alg2 '_iNvals' int2str(iNvals) '_isig' int2str(isig) '_iNtrg' int2str(iNtrg) '_itol' int2str(itol) '.mat'];
itol = 7;
fname2_tol7 = [dir alg2 '_iNvals' int2str(iNvals) '_isig' int2str(isig) '_iNtrg' int2str(iNtrg) '_itol' int2str(itol) '.mat'];
itol = 8;
fname2_tol8 = [dir alg2 '_iNvals' int2str(iNvals) '_isig' int2str(isig) '_iNtrg' int2str(iNtrg) '_itol' int2str(itol) '.mat'];

fname3 = [dir2 'iNvals' int2str(iNvals) '_isig' int2str(isig) '_iNtrg' int2str(iNtrg) '.mat'];
B = load(fname2_tol5);
C = load(fname2_tol7);
D = load(fname2_tol8);

N = Nvals(iNvals);
sig = sigvals(isig);
dim = 2;
ntrgs_per_d = 10^(iNtrg);
xtrgs = equispaced_grid(dim, ntrgs_per_d);

rng(1);
[~,ntrgs] = size(xtrgs);
meastarg = f(xtrgs) + sig*randn(ntrgs,1);

mem_EFGP_5 = (max(B.mginfo.bytes) - min(B.mginfo.bytes))/1024/1024/1024;
mem_EFGP_7 = (max(C.mginfo.bytes) - min(C.mginfo.bytes))/1024/1024/1024;


rms_EFGP_5 = rms(B.ytrgs.mean-meastarg);
rms_EFGP_7 = rms(C.ytrgs.mean-meastarg);
rms_EFGP_8 = rms(D.ytrgs.mean-meastarg);

eepm_new_EFGP_5 = rms(B.ytrgs.mean-D.ytrgs.mean);
eepm_new_EFGP_7 = rms(C.ytrgs.mean-D.ytrgs.mean);

eepm_EFGP_5 = rms(B.y.mean - D.y.mean);
eepm_EFGP_7 = rms(C.y.mean - D.y.mean);

time_EFGP_5 = B.info.cpu_time.total;
time_EFGP_7 = C.info.cpu_time.total;

nxis_EFGP_5 = length(B.info.xis);
nxis_EFGP_7 = length(C.info.xis);

niter_EFGP_5 = B.info.iter;
niter_EFGP_7 = C.info.iter;

mem_FLAM_5 = nan;
mem_FLAM_7 = nan;
eepm_new_FLAM_5 = nan;
eepm_new_FLAM_7 = nan;
time_FLAM_5 = nan;
time_FLAM_7 = nan;
rms_FLAM_5 = nan;
rms_FLAM_7 = nan;
eepm_FLAM_5 = nan;
eepm_FLAM_7 = nan;

if(iNvals <=3)
    A5 = load(fname_tol5);
    A7 = load(fname_tol7);
    mem_FLAM_5 = A5.info.RAM/1024/1024/1024;
    mem_FLAM_7 = A7.info.RAM/1024/1024/1024;

    eepm_new_FLAM_5 = rms(A5.ytrgs.mean-D.ytrgs.mean);
    eepm_new_FLAM_7 = rms(A7.ytrgs.mean-D.ytrgs.mean);

    time_FLAM_5 = A5.info.cpu_time.total;
    time_FLAM_7 = A7.info.cpu_time.total;
    A5.y.mean = meas - A5.alpha*sigmasq;
    A7.y.mean = meas - A7.alpha*sigmasq;
    
    
    
    rms_FLAM_5 = rms(A5.ytrgs.mean-meastarg);
    rms_FLAM_7 = rms(A7.ytrgs.mean-meastarg);
   
    eepm_FLAM_5 = rms(A5.y.mean-D.y.mean);
    eepm_FLAM_7 = rms(A7.y.mean-D.y.mean);

end
save(fname3,'sig','N','nxis_EFGP_5','nxis_EFGP_7','niter_EFGP_5','niter_EFGP_7', ...
       'niter_EFGP_5','niter_EFGP_7','time_EFGP_5','time_EFGP_7','time_FLAM_5', ...
       'time_FLAM_7','mem_EFGP_5','mem_EFGP_7','mem_FLAM_5','mem_FLAM_7',...
       'eepm_EFGP_5','eepm_EFGP_7','eepm_FLAM_5','eepm_FLAM_7',...
       'eepm_new_EFGP_5','eepm_new_EFGP_7','eepm_new_FLAM_5','eepm_new_FLAM_7',...
       'rms_EFGP_5','rms_EFGP_7','rms_FLAM_5','rms_FLAM_7');


