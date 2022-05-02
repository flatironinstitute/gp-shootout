% test 2D rand data with a big hole in it, to test ill-cond, fill-in,
% large-beta instability wrt opts.tol, sigma->0, incons data, etc.
% Currently does a tol-sweep and plots convergence, timing, iters, etc, vs tol.
% Barnett 4/19/22

clear; verb = 1;
dim = 2;        % fixed
N = 1e6;        % data problem size (large, or <=1e4 compares to naive)
M = 1e4;        % test targets covering all incl hole
l = 0.1;        % kernel scale rel to domain L=1. Larger fills hole better.
sigma = .01;    % used to regress. sigma->0 doesn't cause iter # blowup.
sigmadata = 1.0*sigma;   % meas noise: 1.0 for consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified

fprintf('\n2D rand-hole EFGP, l=%.3g, sigma=%.3g...\n',sigma)
rng(1);
hsiz = 0.3;                % hole half-size, must be in (0,1/2),  eg 0.3 std
x = rand(dim,ceil(N/(0.25-hsiz^2))) - 0.5;      % buncha x in [-1/2,1/2]^2
x = x(:,max(abs(x),[],1)>hsiz);      % exclude a central box
x = x(:,1:N);                 % keep just N
xtrg = rand(dim,M) - 0.5;
unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
wavevec = freqdata*unitvec;    % col vec
f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
truemeas = f(x);
truetrg = f(xtrg);
meas = truemeas + sigmadata*randn(N,1);   % noisy data
var0 = var(truemeas);      % use correct prior variance, close to 1 anyway

% choose kernel...
ker = SE_ker(dim,l,var0);         % allows to push to small tol
%nu = 5/2; ker = Matern_ker(dim,nu,l,var0);
tols = 10.^-[2:9];     % conv test: [2:6] for Matern-5/2 or less for smaller nu

nt = numel(tols);
ymeans = nan(N,nt); ytmeans = nan(M,nt);
for i=1:nt, opts.tol = tols(i);       % convergence study..........
  fprintf('EFGP...  tol=%.3g\n',opts.tol);
  %profile clear; profile on
  [y, ytrg, info] = EFGP(x, meas, sigma^2, ker, xtrg, opts);
  %profile off; profile viewer        % yes, it's almost all fftn :)
  ymeans(:,i) = y.mean; ytmeans(:,i) = ytrg.mean;   % save pred vecs
  its(i) = info.iter; tims(i) = info.cpu_time.total; ms(i) = numel(info.xis);
  fprintf('%d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
  fprintf('CPU time (s):'); fprintf('\t%.3g',info.cpu_time.total); fprintf('\n');
  [~,jt] = sort(xtrg(1,:).^2+xtrg(2,:).^2); jt = jt(1);   % targ pt nearest 0
  j=1; fprintf('\t conv at pts y.mean(%d)=%.10g\tytrg.mean(%d)=%.10g\n',j,y.mean(j),jt,ytrg.mean(jt));
  if N<=1e4   % make sure we're computing gp regression accurately...
    fprintf('naive_gp...\n'); % run O(n^3) naive gp regression
    [yn, ytrgn, ~] = naive_gp(x, meas, sigma^2, ker, xtrg, opts);
    fprintf('\trms EFGP vs naive @ train:   %.3g\n', rms(y.mean-yn.mean))
    fprintf('\trms EFGP vs naive @ test:    %.3g\n', rms(ytrg.mean-ytrgn.mean))
  end
end                  % ..............
fprintf('y mean @ data x: rms err vs meas data   %.3g\t(cf sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
fprintf('y mean @ xtrg: rms err vs true f        %.3g\t(expect O(1) for big hole)\n', rms(ytrg.mean-truetrg))

for i=1:nt-1, rs(i) = rms(ymeans(:,i)-ymeans(:,end));  % calc self-conv of vecs
  rts(i) = rms(ytmeans(:,i)-ytmeans(:,end)); end
figure; subplot(2,1,1); loglog(tols(1:end-1), [rs;rts], '+-'); xlabel('tol')
hold on; plot(tols(1:end-1),tols(1:end-1),'r:');
title(sprintf('self-conv: N=%d \\sigma=%.g l=%.g %s',N, sigma, l, ker.fam));
legend('rms y self-conv','rms ytrg self-conv','tol'); axis tight;
subplot(2,1,2); loglog(tols, [its;tims;ms], '+-'); legend('iter counts','CPU time (s)','m modes per dim'); xlabel('tol'); axis tight;

if verb>1, figure; s = 5.0; % blob size.    Spatial plots
  ii = 1:min(3e4,N);   % plot inds
  subplot(2,2,1); scatter(x(1,ii),x(2,ii),s,meas(ii),'filled');
  caxis([-1 1]); c = caxis; axis equal tight; title('y data at x');
  subplot(2,2,2); scatter(x(1,ii),x(2,ii),s,y.mean(ii),'filled');
  caxis(c); axis equal tight; title('GPR mean y at x datapts');
  subplot(2,2,3); scatter(xtrg(1,:),xtrg(2,:),s,truetrg,'filled');
  caxis(c); axis equal tight; title('true f at xtrg');
  subplot(2,2,4); scatter(xtrg(1,:),xtrg(2,:),s,ytrg.mean,'filled');
  caxis(c); axis equal tight; title('GPR mean y at xtrg');
end


% Conclusions:
% hole is not a problem, doesn't cause lots more iters for N<1e4
% But, for large N~1e6, small sigma, find iters grows faster than log(1/eps).
% tol = 1e-5 gives 6e-4 rel l2 err in ytrg.mean (comparing vs tol=1e-6)

% If f too osc (not compatible w/ kernel l: eg freqdata=5 but l=0.2),
% beta blows up and can fail to
% match y meas well, even if tol=1e-6 and only 200 iters, and sigma=0.1.

% Here's a param set causing O(1) changes in ytrg.mean for 

% would |beta| >> 1 indicate max lik (needs det) smaller than optimal?

% iter count is a problem w/ large N, *not* to do with existence/size of hole!

% sigma large leads to wavy fill-in, but good conv @ targs in hole
% misspecified sigma not much problem for stability, |beta| small.


% Self-convergence plots: (SE is quick way to conv-test ill-cond lin sys)

% At sig=1, l=0.1 SE: ytrg rmse gets worse w/ N like tol.N/1e4, for N>1e4
% Whereas, y mean rmse stays about tol to 10tol even for big N.

% at sig=1e-2, l=0.1 SE: even at N=1e4, y rmse ~ (100-1e3)tol, and iff there's
% a hole then ytrg is 1e4.tol.

% at N=1e6 & tol<1e-5, my h=1/(L+Ltime) is not converged - loses 1 digit :(
% Why?

% fixing quadtol and nuffttol at converged vals (eg, 1e-16), then for sig=1e-2
% and N=1e6, the # iters grows as sqrt(1/cgtol), suggesting 1/k^2 iter behavior,
% not geometric c^k CG conv.  This is what causes 1/sqrt(tol) instead of O(tol)
% self-conv. (here l=0.1, std, and it's indep of if hsiz = 0 or 0.3).
