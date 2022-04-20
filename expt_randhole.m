% test 2D rand data with a big hole in it, to test ill-cond, fill-in,
% instability wrt opts.tol, etc.
% Barnett 4/19/22

clear
dim = 2;        % fixed
N = 1e6;        % data problem size (large, or <=1e4 compares to naive)
M = 1e4;        % test targets covering all incl hole
l = 0.1;        % kernel scale rel to domain L=1. Larger fills hole better.
sigma = 0.1;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified

fprintf('\n2D rand-hole EFGP, l=%.3g, sigma=%.3g...\n',sigma)
rng(1);
hsiz = 0.3;                % hole half-size, must be in (0,1/2),  eg 0.25
x = rand(dim,ceil(N*2/(0.25-hsiz^2))) - 0.5;      % x in [-1/2,1/2]^2
x = x(:,max(abs(x),[],1)>hsiz);      % exclude a central box
x = x(:,1:N);                 % keep just N
xtrg = rand(dim,M) - 0.5;
unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
wavevec = freqdata*unitvec;    % col vec
f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
truemeas = f(x);
truetrg = f(xtrg);
meas = truemeas + sigmadata*randn(N,1);   % noisy data
var0 = var(truemeas);      % prior variance, O(1)
% choose kernel...
ker = SE_ker(dim,l,var0);         % allows to push to small tol
%nu = 5/2; ker = Matern_ker(dim,nu,l,var0);

tols = 10.^-[2:6];    % convergence study..........
for tol = tols,  opts.tol = tol;
  fprintf('EFGP...  tol=%.3g\n',opts.tol);
  %profile clear; profile on
  [y, ytrg, info] = EFGP(x, meas, sigma^2, ker, xtrg, opts);
  %profile off; profile viewer        % yes, it's almost all fftn :)
  fprintf('%d iters,\t %d xi-nodes, rms(beta)=%.3g\n',info.iter,numel(info.xis)^dim,rms(info.beta))
  fprintf('CPU times (s):'); fprintf('\t%.3g',info.cpu_time); fprintf('\n');
  % check conv at pts...
  j=1; fprintf('\t conv at pts y.mean(%d)=%.10g\tytrg.mean(%d)=%.10g\n',j,y.mean(j),j,ytrg.mean(j));
  % make sure we're computing gp regression accurately:
  if N<=1e4
    fprintf('naive_gp...\n'); % run O(n^3) naive gp regression
    [yn, ytrgn, ~] = naive_gp(x, meas, sigma^2, ker, xtrg, opts);
    fprintf('        rms efgp vs naive @ train:   %.3g\n', rms(y.mean-yn.mean))
    fprintf('        rms efgp vs naive @ test:    %.3g\n', rms(ytrg.mean-ytrgn.mean))
  end
end                  % ..............
fprintf('y mean @ data x: rms err vs meas data   %.3g\t(cf sigmadata=%.3g)\n', rms(y.mean-meas),sigmadata)
fprintf('y mean @ xtrg: rms err vs true f        %.3g\t(expect O(1) for big hole)\n', rms(ytrg.mean-truetrg))

figure; s = 5.0; % blob size
ii = 1:min(1e5,N);   % plot inds
subplot(2,2,1); scatter(x(1,ii),x(2,ii),s,meas(ii),'filled');
caxis([-1 1]); c = caxis; axis equal tight; title('y data at x');
subplot(2,2,2); scatter(x(1,ii),x(2,ii),s,y.mean(ii),'filled');
caxis(c); axis equal tight; title('GPR mean y at x datapts');
subplot(2,2,3); scatter(xtrg(1,:),xtrg(2,:),s,truetrg,'filled');
caxis(c); axis equal tight; title('true f at xtrg');
subplot(2,2,4); scatter(xtrg(1,:),xtrg(2,:),s,ytrg.mean,'filled');
caxis(c); axis equal tight; title('GPR mean y at xtrg');

% conclusions: hole is not a problem, doesn't cause lots more iters for N<1e4
% But, for large N~1e6, find iters grows faster than log(1/eps).
% tol = 1e-5 gives 6e-4 rel l2 err in ytrg.mean (comparing vs tol=1e-6)

% If f too osc (not compatible w/ kernel l), beta blows up and can fail to
% match y meas well, even if tol=1e-6 and only 200 iters, and sigma=0.1.

% would |beta| >> 1 indicate max lik (needs det) smaller than optimal?
