% Initial experiment to measure error and time on Heaton 2D data.
% The Heaton et al. 2019 paper doesn't state the best/used kernel params for
% this data, but we pulled them out of the github repo.
% Barnett, April 2022.
%
% 'sim' decent sigma~0.6 now (but tol is not indicative, due to ill-cond still):
%  1 sec      tol = 1e-4. (302 its, 30k nodes),   MAE = 0.9
%  6 sec      tol = 3e-5  (420 its, 113k nodes)  MAE = 0.68, ok
%  25 sec for tol=1e-5, (611 its, 380k nodes)   MAE = 0.615,  good.
%  150 sec for tol =3e-6, (814 iters, 1.4M nodes).  MAE = 0.607 (best in paper).
%  (we see y rms resid stil dropping O(1), from 0.42 to 0.35 in last step)

% for 'sat', sigma=0.9 tol=1e-5 gets 0.88 resid, 1.36 MAE, in 22 sec (500 its)

clear; verb = 1;
type = 'sat';
[x,meas,~,xtrg,truetrg] = get_Heatondata(type);
dim = 2;

% kernel params (Heaton et al claim to have fitted) to give as in line 110 of:
%https://github.com/finnlindgren/heatoncomparison/blob/master/Code/FormatData/FormatData.R
mu = 44.49105;      % prior non-zero mean of GP (const in space, apparently)
if strcmp(type,'sim'), sigma = 0.4; else, sigma = 0.9; end  % estim by resid
%sigma = 1.0;        % degree C units; 0.05 way too small, given Wan's MODIS doc
                    % furthermore, sigma=0.05 is unstable (8e3 iters, beta big)
var = 16.40771;     % overall prior variance = covar kernel C(0) (sig rel to it)
l = 1/0.75;   % Matern kernel length in degrees; ignore wrong LON scale :(

nu = 1/2; ker = Matern_ker(dim,nu,l,var); opts.tol = 1e-5;  % tol<1e-5 hurts!
%l = 0.1; ker = SE_ker(dim,l,var); opts.tol = 1e-7;  % try easier kernel

N = numel(meas); Ntrg = numel(truetrg);
fprintf('\ntry mu-shifted EFGP on Heaton: N=%d (Ntrg=%d), sigma=%.3g, tol=%.3g...\n',N,Ntrg,sigma,opts.tol)
tic;   % now do GP regression with non-zero mu handled by subtraction...
[y, ytrg, info] = EFGP(x, meas-mu, sigma^2, ker, xtrg, opts);  % regress (on diff from mu)
y.mean = y.mean + mu; ytrg.mean = ytrg.mean + mu;  % add mu back in
fprintf('done in %.3g sec: \t%d iters, %d xi-nodes, rms(beta)=%.3g\n',toc,info.iter,numel(info.xis)^dim,rms(info.beta))
fprintf('CPU time, total (s)  :'); fprintf('\t%.3g',info.cpu_time.total); fprintf('\n');
fprintf('CPU time, precomp (s):'); fprintf('\t%.3g',info.cpu_time.precomp); fprintf('\n');
fprintf('CPU time, cg (s)     :'); fprintf('\t%.3g',info.cpu_time.cg); fprintf('\n');
fprintf('CPU time, mean (s)   :'); fprintf('\t%.3g',info.cpu_time.mean); fprintf('\n');
% extract a good choice for sigma, around 0.4-0.6... *not* 0.05!
fprintf('rms resid @ x train pts:\t%.3g (degrees C, suggests sigma)\n', rms(y.mean-meas))

MAE = mean(abs(ytrg.mean-truetrg));
fprintf('MAE @ test targs:       \t%.3g (degrees C; sim aims for 0.6)\n', MAE)

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
