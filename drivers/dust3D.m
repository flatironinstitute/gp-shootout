% driver script for an experiment. (copied from self-test of EFGP)
% This does not do the dust line-integral measurement type!
% (which needs quadrature, etc, and breaks the Toeplitz property)
% Barnett 5/22/23.
clear; verb = 1;

N = 1e7;        % problem size
l = 0.05;        % kernel scale (recall domain is cube length 1)
sigma = 0.03;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 5.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
opts.tol = 1e-4;
dim = 3;
unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
wavevec = freqdata*unitvec;    % col vec
f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
[x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);
ker = SE_ker(dim,l);
%ker = Matern_ker(dim,3/2,l);   % harder for EFGP

fprintf('\ndust3D N=%d, sigma=%.3g, tol=%.3g...\n',N,sigma,opts.tol)
[y, ~, info] = EFGP(x, meas, sigma^2, ker, [], opts);     % regress
fprintf('EFGP %d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
fprintf('CPU time (s):\n'); disp(info.cpu_time)
fprintf('y.mean: rms err vs meas data   %.3g\t(should be about sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
% estim ability to average away noise via # pts in the rough kernel support...
fprintf('        rms truemeas pred err  %.3g\t(should be sqrt(l^d.N) better ~ %.2g)\n', rms(y.mean-truemeas),sigmadata/sqrt(l^dim*N))

if verb, figure; siz = 1.0; Np = min(N,1e5);  % max pts to plot
  subplot(1,2,1); scatter3(x(1,1:Np),x(2,1:Np),x(3,1:Np),siz*ones(Np,1),meas(1:Np),'filled');
  caxis([-1 1]); axis equal tight vis3d
  subplot(1,2,2); scatter3(x(1,1:Np),x(2,1:Np),x(3,1:Np),siz*ones(Np,1),y.mean(1:Np),'filled');
  caxis([-1 1]); axis equal tight vis3d
  title(sprintf('dust3D: %dd N=%d',dim,N));
end
