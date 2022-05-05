% Initial experiment to measure GPR error and time on Heaton 2D data.
% The model has non-zero mean (mu) but this is handled by subtracting before
% GPR, then adding back in.
% The Heaton et al. 2019 paper doesn't state the best/used kernel params for
% this data, but we pulled them out of the github repo.
% Barnett, April 2022. Various methods, 5/2/22. RLCM 5/5/22.
%
% old EFGP results.... v sensitive to tol :(
% 'sim' decent sigma~0.6 now (but tol is not indicative, due to ill-cond still):
%  1 sec      tol = 1e-4. (302 its, 30k nodes),   MAE = 0.9
%  6 sec      tol = 3e-5  (420 its, 113k nodes)  MAE = 0.68, ok
%  25 sec for tol=1e-5, (611 its, 380k nodes)   MAE = 0.615,  good.
%  150 sec for tol =3e-6, (814 iters, 1.4M nodes).  MAE = 0.607 (best in paper).
%  (we see y rms resid stil dropping O(1), from 0.42 to 0.35 in last step)
% We see that sigma = 0.4-0.6 is a realistic choice for meas error, *not* 0.05!
% for 'sat', sigma=0.9 tol=1e-5 gets 0.88 resid, 1.36 MAE, in 22 sec (500 its)

clear; verb = 1;  % 0 no figs, 1 figs.
type = 'sat';     % choose: 'sim' (simulated GP) or 'sat' (real data)
[x,meas,~,xtrg,truetrg] = get_Heatondata(type);
dim = 2;  % fix

% kernel params (Heaton et al claim to have fitted) to give as in line 110 of:
%https://github.com/finnlindgren/heatoncomparison/blob/master/Code/FormatData/FormatData.R
mu = 44.49105;      % prior non-zero mean of GP (const in space, apparently)
if strcmp(type,'sim'), sigma = 0.4; else, sigma = 0.9; end  % estim by resid
%sigma = 1.0;        % degree C units; 0.05 way too small, given Wan's MODIS doc
                    % furthermore, sigma=0.05 is unstable (8e3 iters, beta big)
var = 16.40771;     % overall prior variance = covar kernel C(0) (sig rel to it)
l = 1/0.75;   % Matern kernel length in degrees; ignore wrong LON scale :(

nu = 1/2; ker = Matern_ker(dim,nu,l,var);
%l = 0.1; ker = SE_ker(dim,l,var);  % try easier wrong kernel

N = numel(meas); Ntrg = numel(truetrg);
for meths={'EFGP','FLAMGP','SKI','RLCM'}           % -- METHODS LOOP -----
  meth=meths{1};      % make not a 1-element cell array
  fprintf('\nHeaton, method %s: \tN=%d (Ntrg=%d), sigma=%.3g...\n',meth,N,Ntrg,sigma)
  % now do GP regression with non-zero mu handled by subtraction,
  % ie regress on diff from mu....
  switch meth
    case 'EFGP'
      opts = []; opts.tol = 1e-5;  % tol<1e-5 hurts!
      [y, ytrg, info] = EFGP(x, meas-mu, sigma^2, ker, xtrg, opts);
      fprintf('\t%d iters, %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
    case 'FLAMGP'
      opts = []; opts.tol = 1e-5;
      [y, ytrg, info] = FLAMGP(x, meas-mu, sigma^2, ker, xtrg, opts);
    case 'SKI'
      opts = []; opts.grid_size = 300;    % who knows?
      [y, ytrg, info] = SKI(x, meas-mu, sigma^2, ker, xtrg, opts);
    case 'RLCM'
      opts = []; opts.rank = 250;   % higher than 200 needed for decent 2d acc
      [y, ytrg, info] = RLCM(x, meas-mu, sigma^2, ker, xtrg, opts);
  end
  y.mean = y.mean + mu; ytrg.mean = ytrg.mean + mu;  % add mu back in
  fprintf('CPU time (s):\n'); disp(info.cpu_time);
  fprintf('rms resid @ x train pts:\t%.3g (degrees C)\n', rms(y.mean-meas))
  MAE = mean(abs(ytrg.mean-truetrg));
  fprintf('MAE @ test targs:       \t%.3g (degrees C; sim aims for 0.6)\n', MAE)
  
  if verb, figure;
    xa = [x, xtrg];        % gather all nodes
    ya = [meas; truetrg];  % all true vals
    subplot(2,1,1); scatter(xa(1,:),xa(2,:),1,ya); axis equal tight; colorbar
    c = caxis;
    title(sprintf('Heaton: all (train+test) data'));
    subplot(2,1,2); scatter(xtrg(1,:),xtrg(2,:),1,ytrg.mean); axis equal tight;
    caxis(c); colorbar
    title(sprintf('Heaton %s: pred mean targ data',meth)); drawnow
  end
end                                         % --------------------------------


% RESULTS TABLE

% results for 'sat' dataset @ sigma=0.9:
% EFGP @ tol=1e-5:  22sec,  0.877 rms resid,   1.36 MAE test
% FLAMGP tol=1e-5:  26sec,  0.669 rms (v good) 1.37 MAE test
% SKI    grid=300   14sec,  1.32   "           1.86 "   (v inacc; see rect errs)
% RLCM   rank=250   52s  ,  0.673     (v good) 1.42 "         (ok)

% results for 'sim' dataset @ sigma=0.4:
% EFGP @ tol=1e-5:  34sec,  0.414 rms resid,  0.624 MAE test  (good: aim 0.6)
% FLAMGP tol=1e-5:  26sec,  0.226 rms (good), 0.602 MAE test  (good)
% SKI    grid=300   18sec,  0.564   "         0.741 "    (inacc; see rect errs)
% RLCM   rank=250   51s  ,  0.225     (good)  0.634 "         (good)
