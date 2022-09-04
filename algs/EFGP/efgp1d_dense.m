function [beta, xis, ytrg, time_info, A, X, ws] = efgp1d_dense(x, y, sigmasq, ker, eps, xsol, opts)
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
  
    if(nargin == 6), opts = []; end
    k = ker.k; khat = ker.khat;  % get functions, new kernel format
    N = numel(y);
    
    tic_precomp = tic;
    x0 = min([x; xsol]); x1 = max([x; xsol]);
    L = x1-x0;                   % approx domain length
    [xis, h, m] = get_xis(ker, eps, L);

    % khat & quadr weight scaling of Fourier basis functions
    ws = sqrt(khat(xis)' * h);
    t_precomp = toc(tic_precomp);
    
    tic_solve = tic; 
    X = exp(1i * x * xis * 2 * pi);
    A = diag(ws) * (X'* X) * diag(ws);
    beta = (A + sigmasq * eye(m)) \ (diag(ws) * X' * y);
    t_solve = toc(tic_solve);
    
    % tabulate posterior mean
    tic_post = tic;
    Xsol = exp(1i * xsol * xis * 2 * pi) * diag(ws);
    yhat = Xsol * beta;
    ytrg.mean = real(yhat);
    t_post = toc(tic_post);
    
    % if specified by user, compute posterior variance at target points via
    % the the diagonal of F C^{-1} F' where C is the posterior covariance
    % matrix and F is the matrix of basis functions (complex exponentials)
    % tabulated at target points
    do_var = false;
    if(isfield(opts,'get_var'))
        do_var = opts.get_var;
    end
    if do_var
        nsol = numel(xsol);
        ytrg.var = zeros(nsol, 1);
        c = (A / sigmasq) + eye(m);
        c_inv = inv(c);
        fs = exp(1i * 2 * pi * xsol * xis) * diag(ws);
        ytrg.var = real(diag(fs * c_inv * fs'));
    end

    time_info = [t_precomp, t_solve, t_post, t_post - t_precomp];
end
