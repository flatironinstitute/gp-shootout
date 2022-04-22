function [beta, xis, yhat, time_info, A, ws] = function_space1d_dense(x, y, sigmasq, ker, eps, xsol)
% FUNCTION_SPACE1D_DENSE   equispaced Fourier GP regression in 1D by
% constructing the matrix of the linear system and doing a dense solve. 
%
% [beta, xis, yhat, iter, time_info] = function_space1d_dense(x, y, sigmasq, ker, eps, xsol) 
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
% xsol   - locations at which to evaluate posterior mean
%
% Outputs:
% beta - vector of Fourier basis weights (not really for the user)
% xis  - Fourier freqs used (not really for the user)
% yhat - posterior means at xsol ordinates   <- the only user output
% time_info  - diagnostic list of timings
% X - matrix of linear system
% ws - scaling of each complex exponential
%
% To test this routine see: EFGP
  
    k = ker.k; khat = ker.khat;  % get functions, new kernel format
    N = numel(y);
    
    tic_precomp = tic;
    x0 = min(x); x1 = max(x);
    L = x1-x0;                   % approx domain length *** could check xtrg too?
    [xis, h, m] = get_xis(ker, eps, L);
    % center all coords for NUFFTs domain, then do 2pi.h ("tph") rescaling...
    xcen = (x1+x0)/2;
    tphx = 2*pi*h*(x - xcen);
    tphxsol = 2*pi*h*(xsol - xcen);
    
    % khat & quadr weight scaling of Fourier basis functions
    ws = sqrt(khat(xis)' * h);
    t_precomp = toc(tic_precomp);
    
    tic_solve = tic; 
    X = exp(1i * x * xis * 2 * pi);
    A = diag(ws) * (X'* X) * diag(ws);
    beta = (A + sigmasq * eye(m)) \ (diag(ws) * X' * y);
    t_solve = toc(tic_solve);
    
    % tabulate solution using fft
    tic_post = tic;
    tmpvec = ws .* beta;
    nuffttol = eps / 10;   % nufft is fast, so keep its errors insignificant
    yhat = finufft1d2(tphxsol, +1, nuffttol, tmpvec);
    
    % tabulate posterior mean
    Xsol = exp(1i * xsol * xis * 2 * pi) * diag(ws);
    yhat = Xsol * beta;
    yhat = real(yhat);
    t_post = toc(tic_post);
    
    time_info = [t_precomp, t_solve, t_post, t_post - t_precomp];
end
