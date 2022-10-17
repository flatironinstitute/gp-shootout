rng(1);
sigmatrue = 0.1;

% set directory for saving data
dir = "~/ceph/gp-shootout/big_example/data";

% 2d
N = 1e9;
dim = 2;
wavevec = [4; 3];
phase = 1.3;
f = @(x) cos(2*pi*x'* wavevec + phase);
x = rand(dim,N);
truemeas = f(x);

sigmas = [0.1; 0.1*sqrt(1e7/3e6); 0.1*sqrt(3e7/3e6); 0.1*sqrt(1e8/3e6); 0.1*sqrt(1e9/3e6)];
nsig = length(sigmas);
for i=1:nsig
  disp(i);
  meas = truemeas(:) + sigmas(i)*randn(N,1);
  fstr = ['meas_2d_isig' int2str(i) '.bin'];
  disp(fstr);
  fid = fopen(fullfile(dir,fstr),'w');
  fwrite(fid,meas,'double');
  fclose(fid);
end

tic,
fid = fopen(fullfile(dir,'x_2d.bin'),'w');
fwrite(fid,x,'double');
toc;
fclose(fid);

