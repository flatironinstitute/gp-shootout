% Initial experiment to measure error and time on Heaton 2D data.
% Sadly, the Heaton et al. 2019 paper doesn't give the best kernel params for
% this data. Also we have a major mean=0 problem!
% Barnett 4/11/22

clear; verb = 1;

load heaton19sim   % path gets it: goes to workspace: x meas xtrg truetrg
N = numel(meas);
Ntrg = numel(truetrg);
l = 0.05;        % SE kernel scale - who knows?
sigma = 0.1;    % used to regress - should be meas err - at least 0.5 in MODIS
opts.tol = 1e-3;   % crude
dim = 2;
%ker = SE_ker(dim,l);  % or
ker = Matern_ker(dim,1/2,l);

fprintf('\ntry EFGP on Heaton: N=%d (Ntrg=%d), sigma=%.3g, tol=%.3g...\n',N,Ntrg,sigma,opts.tol)
tic;
[y, ytrg, info] = EFGP(x, meas, sigma^2, ker, xtrg, opts);     % regress
y.mean = real(y.mean); ytrg.mean = real(ytrg.mean);
fprintf('done in %.3g sec: \t%d iters,\t %d xi-nodes, rms(beta)=%.3g\n',toc,info.iter,numel(info.xis)^dim,rms(info.beta))
fprintf('CPU time breakdown (s):'); fprintf('\t%.3g',info.cputime); fprintf('\n');

MAE = mean(abs(ytrg.mean-truetrg));
fprintf('MAE @ test targs: %.3g (degrees Celsius or K; aiming for 0.6)\n', MAE)
% well, since mean = 0, it's a terrible predictor!!   *** to investigate
% Q: What is the right GP regression test here - to fit kernel params?

if verb, figure;
  xa = [x, xtrg];        % gather all nodes
  ya = [meas; truetrg];  % all true vals
  subplot(1,2,1); scatter(xa(1,:),xa(2,:),1,ya); axis equal tight; colorbar
  c = caxis;
  title('all (train+test) data');
  ypred = [meas; ytrg.mean];  % train + pred mean at ytrg
  subplot(1,2,2); scatter(xa(1,:),xa(2,:),1,ypred); axis equal tight;
  caxis(c); colorbar
  title('train + pred mean targ data');
end
