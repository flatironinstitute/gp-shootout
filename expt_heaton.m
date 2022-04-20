% Initial experiment to measure error and time on Heaton 2D data.
% The Heaton et al. 2019 paper doesn't state the best/used kernel params for
% this data, but we pulled them out of the github repo.
% Barnett, April 2022.

clear; verb = 1;
[x,meas,~,xtrg,truetrg] = get_Heatondata('sim');

% kernel params (Heaton claim to have fitted) to give as in line 110 of:
%https://github.com/finnlindgren/heatoncomparison/blob/master/Code/FormatData/FormatData.R
mu = 44.49105;      % prior non-zero mean of GP (const in space, by why not try linear, since obvious SW-increasing trend?)
sigma = 1;       % in degree C units; 0.05 too small given Wan's MODIS doc
var = 16.40771;     % overall prior variance = covar kernel C(0)
l = 1/0.75;   % kernel lengthscale in "degree" units; ignore wrong LON scale :(
dim = 2; nu = 1/2;
%ker = Matern_ker(dim,nu,l,var);
ker = SE_ker(dim,l,var);
opts.tol = 1e-7;    % alg param: has big effect on CPU time for nu=1/2

N = numel(meas); Ntrg = numel(truetrg);
fprintf('\ntry mu-shifted EFGP on Heaton: N=%d (Ntrg=%d), sigma=%.3g, tol=%.3g...\n',N,Ntrg,sigma,opts.tol)
tic;   % now do GP regression with non-zero mu handled by subtraction...
[y, ytrg, info] = EFGP(x, meas-mu, sigma^2, ker, xtrg, opts);  % regress (on diff from mu)
y.mean = y.mean + mu; ytrg.mean = ytrg.mean + mu;  % add mu back in
fprintf('done in %.3g sec: \t%d iters, %d xi-nodes, rms(beta)=%.3g\n',toc,info.iter,numel(info.xis)^dim,rms(info.beta))
fprintf('CPU time breakdown (s):'); fprintf('\t%.3g',info.cpu_time); fprintf('\n');

MAE = mean(abs(ytrg.mean-truetrg));
fprintf('MAE @ test targs: %.3g (degrees Celsius or K; aiming for 0.6)\n', MAE)

if verb, figure;
  xa = [x, xtrg];        % gather all nodes
  ya = [meas; truetrg];  % all true vals
  subplot(1,2,1); scatter(xa(1,:),xa(2,:),1,ya); axis equal tight; colorbar
  c = caxis;
  title('all (train+test) data');
  subplot(1,2,2); scatter(xtrg(1,:),xtrg(2,:),1,ytrg.mean); axis equal tight;
  caxis(c); colorbar
  title('pred mean targ data');
end

% *** to do, check tol effect on y.mean, to check convergence
