function [beta, xis, yhat, iter, time_info] = function_space1d(x, y, sigmasq, ker, eps, xsol)
% FUNCTION_SPACE1D   fast equispaced Fourier NUFFT-based GP regression in 1D
%
% [beta, xis, yhat, iter, time_info] = function_space1d(x, y, sigmasq, ker, eps, xsol) 
% performs Gaussian process regression in 1d using equispaced
% Fourier representations of Gaussian processes and fast algorithms for
% performing regression. 
%
% Inputs:
% x      - N x 1 array of location of observations
% y      - N x 1 array of (noisy) observations
% sigmasq - residual variance for GP regression
% ker    - struct with ker.k is the covariance kernel and ker.khat is its
% Fourier transform
% eps    - truncate covariance kernel in time and Fourier domains when values
% of functions reach eps
% xsol   - locations at which to evaluate posterior mean
%
% Outputs:
% beta - vector of Fourier basis weights (not really for the user)
% xis  - Fourier freqs used (not really for the user)
% yhat - posterior means at xsol ordinates   <- the only user output
% iter - diagnostics from CG
% ts   - diagnostic list of timings
%
% To test this routine see: EFGP
  
  k = ker.k; khat = ker.khat;    % get functions, new kernel format
  N = numel(y);

    tic_precomp = tic;
    tmax = 1;
    dim = 1;
    xis = get_xis(dim, ker, eps, tmax);
    h = xis(2) - xis(1);
    m = numel(xis);

    % set scaling of basis functions
    ws = sqrt(khat(xis)' * h);
    
    % construct first row and column of toeplitz matrix for fast apply
    tol = eps / 10; % nufft is fast, so just make sure we don't incur errors
    c = ones(N, 1);
    XtXrow = finufft1d1(x*2*pi*h, c, +1, tol, 2*m-1)'; 
    Gf = fftn(XtXrow.');
    
    % construct rhs with fft
    isign = -1;
    tol = eps / 10; % nufft is fast, so just make sure we don't incur errors
    rhs = finufft1d1(2*pi*x*h, y, isign, tol, m);
    rhs = ws .* rhs;
    
    t_precomp = toc(tic_precomp);
    
    % solve linear system with conjugate gradient
    Afun = @(a) ws .* Afun2(Gf, ws .* a) + sigmasq .* a; 
    
    tic_cg = tic; 
    [beta,flag,relres,iter,resvec] = pcg(Afun, rhs, eps, m);
    t_cg = toc(tic_cg);
    
    % tabulate solution using fft
    tol = eps / 10; % nufft is fast, so just make sure we don't incur errors
    tmpvec = ws .* beta;
    tic_post = tic;
    yhat = finufft1d2(2*pi*h*xsol, +1, tol, tmpvec);
    t_post = toc(tic_post);

    time_info = [t_precomp, t_cg, t_post];

    % convert to real
    yhat = real(yhat);
end


function [v2] = Afun2(Gf, a)
% this function is used for performaing a fast matrix multiply 
% in conjugate gradient. it takes an fft of a and multiplies 
% in frequency domain by Gf and then converts back to time/spatial domain
    m = numel(a);
    af = fftn(a, size(Gf));
    vft = af .* Gf;
    vft = ifftn(vft);
    v2 = vft(m:end);
end