rng(1);
l = 0.1;          % SE kernel scale
sigmatrue = 0.3;  % used to regress
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified
for dim = 1:3
    fprintf('\n%g dimension(s) \n\n', dim);
    for N = [1e2, 1e3, 1e4]
        fprintf('dim=%g...\n', dim)
        unitvec = randn(dim, 1); unitvec = unitvec/norm(unitvec);
        wavevec = freqdata * unitvec;    % col vec
        f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
        [x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);
        ker = SE_ker(dim,l);
        sigmasq = sigmatrue^2;
        ntrgs_per_d = 100;
        xtrgs = equispaced_grid(dim, ntrgs_per_d);
        
        opts.tol = 1e-4;
        [y, ytrgs, info] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
        
        opts.tol = opts.tol / 50;
        [ytrue, ytrgs_true, info_true] = EFGP(x, meas, sigmasq, ker, xtrgs, opts);
        fprintf('EFGP rms at targets %.3g, time: %.3g\n', rms(ytrgs.mean-ytrgs_true.mean), info.cpu_time.total);
        fprintf('rms  at targets %.3g\n', rms(ytrgs.mean-ytrgs_true.mean));
        fprintf('linf at targets %.3g\n', max(abs(ytrgs.mean-ytrgs_true.mean)));
        fprintf('total time   %.3g\n', info.cpu_time.total);
        fprintf('precomp time %.3g\n', info.cpu_time.precomp);
        fprintf('cg time      %.3g\n', info.cpu_time.cg);
        fprintf('mean time    %.3g\n', info.cpu_time.mean);
        fprintf('mean/target  %.3g\n', info.cpu_time.mean / (ntrgs_per_d^dim));
        
        %%%fprintf("$ 10^{%d}$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $ %d $ & $ %d $ \\\\ \n", ...
    end
end
