function [beta, xis, ytrg, iter, time_info] = efgp1d(x, y, sigmasq, ker, eps, xtrgs, do_var)
% EFGP1D   fast equispaced Fourier NUFFT-based GP regression in 1D
%
% [beta, xis, ytrg, iter, time_info] = efgp1d(x, y, sigmasq, ker, eps, xtrgs) 
% performs Gaussian process regression in 1d using equispaced
% Fourier representations of Gaussian processes and fast algorithms for
% performing regression. 
%
% Inputs:
% x      - N x 1 array of location of observations
% y      - N x 1 array of (noisy) observations
% sigmasq - residual variance for GP regression
% ker    - struct with ker.k is the covariance kernel and ker.khat is its
%          Fourier transform
% eps    - truncate covariance kernel in time and Fourier domains when values
%          of functions reach eps
% xtrgs   - locations at which to evaluate posterior mean
%
% Outputs:
% beta - vector of Fourier basis weights (not really for the user)
% xis  - Fourier freqs used (not really for the user)
% ytrg - ytrg.mean - posterior means at xtrgs ordinates   <- the only user output
%        ytrg.var - [optional] posterior variances at xtrgs ordinates
% iter - diagnostics from CG
% time_info  - diagnostic list of timings
%
% To test this routine see: EFGP
    
    k = ker.k; khat = ker.khat;  % get functions, new kernel format
    N = numel(y);
    
    tic_precomp = tic;
    x0 = min([x; xtrgs]); x1 = max([x; xtrgs]);
    L = x1-x0;                   % approx domain length *** could check xtrg too?
    [xis, h, m] = get_xis(ker, eps, L);
    % center all coords for NUFFTs domain, then do 2pi.h ("tph") rescaling...
    xcen = (x1+x0)/2;
    tphx = 2*pi*h*(x - xcen);
    tphxtrgs = 2*pi*h*(xtrgs - xcen);
    
    % khat & quadr weight scaling of Fourier basis functions
    ws = sqrt(khat(xis)' * h);
    
    % construct first row and column of toeplitz matrix for fast apply
    nuffttol = eps / 10;   % nufft is fast, so keep its errors insignificant
    c = complex(ones(N, 1));      % unit weights
    XtXcol = finufft1d1(tphx, c, -1, nuffttol, 2*m-1);
    Gf = fftn(XtXcol);          % no conj; equiv to isign=-1
    
    % construct rhs = X^*y, with NUFFT
    isign = -1;
    rhs = finufft1d1(tphx, y, isign, nuffttol, m);
    rhs = ws .* rhs;                         % col vecs
    t_precomp = toc(tic_precomp);
    
    % solve linear system (X^*X + sigma^2)beta = rhs with conjugate gradient
    Afun = @(a) ws .* Afun2(Gf, ws .* a) + sigmasq .* a; 
    
    tic_cg = tic; 
    [beta,flag,relres,iter,resvec] = pcg(Afun, rhs, eps, 3*m);  % solve beta vec
    t_cg = toc(tic_cg);
    
    % tabulate solution using fft
    tmpvec = ws .* beta;
    tic_post = tic;
    yhat = finufft1d2(tphxtrgs, +1, nuffttol, tmpvec);
    ytrg.mean = real(yhat);
    t_post = toc(tic_post);
    
    % if specified by user, compute posterior variance at target points. at
    % each target point evaluate the quadratic form f' C^{-1} f where C is
    % the posterior covariance and f is a vector of basis functions
    % (complex exponentials) tabulated at one target point
    if do_var
        eps_var = max(1e-3, eps); % don't need high precision
        Afun_var = @(a) ws .* Afun2(Gf, ws .* a)/sigmasq + a; 
        ntrgs = numel(xtrgs);
        ytrg.var = zeros(ntrgs, 1);
        for i=1:ntrgs
            rhs = exp(1i * tphxtrgs(i) * (-(m-1)/2:1:(m-1)/2)) * diag(ws);
            [beta,flag,relres,iter,resvec] = pcg(Afun_var, rhs', eps_var, 3*m);
            ytrg.var(i) = real(rhs * beta);
        end
    end

    % timings
    time_info = [t_precomp, t_cg, t_post];
    time_info = [time_info, sum(time_info)];

end


function [v2] = Afun2(Gf, a)
% this function performs a fast multiply by Toeplitz matrix, needed for
% conjugate gradient. it takes an fft of a and multiplies 
% in "coeff frequency" domain by Gf and then converts back to "coeff" domain.
    m = numel(a);
    af = fftn(a, size(Gf));
    vft = af .* Gf;
    vft = ifftn(vft);
    v2 = vft(m:end);
end
