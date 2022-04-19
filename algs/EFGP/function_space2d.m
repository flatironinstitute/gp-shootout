function [beta, xis, yhat, iter, time_info] = function_space2d(xs, y, sigmasq, ker, eps, xsols)
% FUNCTION_SPACE2D   equispaced Fourier NUFFT-based GP regression in 2D
%
% [beta, xis, yhat, iter, time_info] = function_space2d(xs, y, sigmasq, ker, eps, xsols)
% performs Gaussian process regression in 2d using equispaced
% Fourier representations of Gaussian processes and fast algorithms for
% performing regression. 
%
% Inputs:
% xs - N x 2 array of location of observations
% y - N x 1 array of (noisy) observations
% sigmasq - residual variance for GP regression
% ker - struct with ker.k is the covariance kernel and ker.khat is its
% Fourier transform
% eps - truncate covariance kernel in time and Fourier domains when values
% of functions reach eps
% xsol - locations at which to evaluate posterior mean
%
% Outputs:
% beta - vector of Fourier basis weights (not really for the user)
% xis  - Fourier freqs used (not really for the user)
% yhat - posterior means at xsol ordinates   <- the only user output
% iter - diagnostics from CG
% time_info   - diagnostic list of timings
%
% For test see EFGP 

% get kernel functions
  k = ker.k; khat = ker.khat;

    % support of functionin time domain
    tic_precomp = tic;
    tmax = 1;
    xis = get_xis(ker, eps, tmax);
    h = xis(2) - xis(1);
    m = numel(xis);
    [xis_xx, xis_yy] = meshgrid(xis);

    const = h;
    rs = sqrt(xis_xx.^2 + xis_yy.^2);
    ws = sqrt(khat(rs)) * const;
    
    % precomputation for fast apply of X*X
    Gf = getGf(eps, xs, xis);
    
    % conjugate gradient
    ws_flat = reshape(ws, m^2, 1);
    Afun = @(a) ws_flat .* apply_xtx(Gf, ws_flat .* a) + sigmasq .* a;
    
    isign = -1;
    tol = eps / 10;
    rhs = finufft2d1(2*pi*xs(:,1)*h, 2*pi*xs(:,2)*h, y, isign, tol, m, m);
    rhs = reshape(rhs.' .* ws, [], 1); 
    
    t_precomp = toc(tic_precomp);
    % solve linear system
    tic_cg = tic;
    [beta,flag,relres,iter,resvec] = pcg(Afun,rhs,eps,m^2);
    t_cg = toc(tic_cg);
    
    % evaluate posterior mean 
    tic_post = tic;
    tmpvec = ws .* reshape(beta, [m, m]);
    tmpvec = tmpvec.';
    isign = +1;
    tol = eps / 10;
    yhat = finufft2d2(2*pi*h * xsols(:,1), 2*pi*h * xsols(:,2),isign,tol,tmpvec);
    t_post = toc(tic_post);

    % package timings
    time_info = [t_precomp, t_cg, t_post];
    time_info = [time_info, sum(time_info)];

    % convert to real
    yhat = real(yhat);

end


function [Gf] = getGf(eps, xs, z)
    N = length(xs);
    m = length(z);
    % determine all the differences where kernel is to be evaluated...
    % ...for both x and y directions 
    h = z(2) - z(1);

    % parameters for fft
    c = 0i + ones(N, 1);
    isign = -1;
    tol = eps/10;
    
    o.modeord = 1;
    Gf = finufft2d1(2*pi*xs(:,1)*h,2*pi*xs(:,2)*h,c,isign,tol,2*m-1, 2*m-1, o);
    Gf = fftn(Gf.');
end


function [v] = apply_xtx(Gf, b)
    n = sqrt(length(b));
    b = reshape(b, n, n);
    vft = fftn(b,size(Gf));
    vft = vft.*Gf;
    vft = ifftn(vft);
    v = vft(1:n,1:n);
    v = reshape(v,[],1);
end

