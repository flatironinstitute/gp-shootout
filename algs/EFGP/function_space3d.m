function [beta, xis, yhat, iter, time_info] = function_space3d(x, y, sigmasq, ker, eps, xsol)
% FUNCTION_SPACE3D   equispaced Fourier NUFFT-based GP regression in 3D
%
% [beta, xis, yhat, iter, time_info] = function_space2d(xs, y, sigmasq, ker, eps, xsols)
% performs Gaussian process regression in 3D using equispaced
% Fourier representations of Gaussian processes and fast algorithms for
% performing regression.
%
% Inputs:
% x      - N x 3 array of location of observations
% y      - N x 1 array of (noisy) observations
% sigmasq - residual variance for GP regression
% ker - struct with ker.k is the covariance kernel and ker.khat is its
%       Fourier transform
% eps - truncate covariance kernel in time and Fourier domains when values
%         of functions reach eps
% xsol - n*3 location coords at which to evaluate posterior mean
%
% Outputs:
% beta - vector of Fourier basis weights (not really for the user)
% xis  - 1D grid of Fourier freqs used (not really for the user)
% yhat - posterior means at xsol ordinates   <- the only user output
% iter - diagnostics from CG
% time_info   - diagnostic list of timings
%
% For test see EFGP 

% Note: no attempt to exploit cuboid box done; enclosing cube used for L

  tic_precomp = tic;
  x0 = min(x); x1 = max(x);   % both row 2-vectors
  L = max(x1-x0);    % worst-axis domain length *** could check xsol too?
  [xis h m] = get_xis(ker, eps, L);
  [xis_xx, xis_yy, xis_zz] = ndgrid(xis,xis,xis);        % assumes isotropic
  % center all coords for NUFFTs domain, then do 2pi.h ("tph") rescaling...
  xcen = (x1+x0)/2;                    % row vec
  tphx = 2*pi*h*(x - xcen);            % note broadcast over rows
  tphxsol = 2*pi*h*(xsol - xcen);      % "

  % weights of Fourier basis funcs
  rs = sqrt(xis_xx.^2 + xis_yy.^2 + xis_zz.^2);
  dim = 3; ws = sqrt(khat(rs) * h^dim);
    
  % precomputation for fast apply of X*X
  nuffttol = eps/10;
  Gf = getGf3d(nuffttol, tphx, m);
  
    % conjugate gradient
    ws_flat = reshape(ws, m^3, 1);
    Afun = @(a) ws_flat .* apply_xtx3d(Gf, ws_flat .* a) + sigmasq .* a;
    
    isign = -1;
    a = 2*pi*xs(:,2)*h;
    b = 2*pi*xs(:,1)*h;
    c = 2*pi*xs(:,3)*h;
    tol = eps/10;
    rhs2 = finufft3d1(a, b, c, y, isign, nuffttol, m, m, m);
    rhs2 = reshape(rhs2 .* ws, [], 1); 
    t_precomp = toc(tic_precomp);
    
    % solve linear system
    tic_cg = tic; 
    tol = eps/10;
    [beta,flag,relres,iter,resvec] = pcg(Afun, rhs2, tol, m^3);
    t_cg = toc(tic_cg);

    % evaluate posterior mean 
    tic_post = tic;
    tmpvec = ws .* reshape(beta, [m, m, m]);
    isign = +1;
    tol = eps/10;
    a = 2*pi*h * xsols(:,2);
    b = 2*pi*h * xsols(:,1);
    c = 2*pi*h * xsols(:,3);
    yhat = finufft3d2(a, b, c, isign, nuffttol, tmpvec);
    t_post = toc(tic_post);

    % package times for output
    time_info = [t_precomp, t_cg, t_post];
    time_info = [time_info, sum(time_info)];

    % convert to real
    yhat = real(yhat);

end


function [Gf] = getGf3d(nuffttol, xs, z)
    N = length(xs);
    m = length(z);
    % determine all the differences where kernel is to be evaluated...
    % ...for both x and y directions 
    xs_tmp = z - z(1);
    xs_lrg = [xs_tmp,-xs_tmp(end:-1:2)];
    
    % precompute elements of convolution vector for quick apply
    [ds_xx, ds_yy, ds_zz] = ndgrid(xs_lrg);
    ds_xx_flat = reshape(ds_xx, [], 1);
    ds_yy_flat = reshape(ds_yy, [], 1);
    ds_zz_flat = reshape(ds_zz, [], 1);

    % parameters for fft
    c = complex(ones(N, 1));            % unit strengths
    isign = -1;
    x = 2*pi*xs(:,1);
    y = 2*pi*xs(:,2);
    z = 2*pi*xs(:,3);
    Gf = finufft3d3(x,y,z,c,isign,nuffttol,ds_xx_flat,ds_yy_flat,ds_zz_flat);
    Gf = reshape(Gf, m*2 - 1, m*2 - 1, m*2 - 1);
    Gf = fftn(Gf);
end


function [v] = apply_xtx3d(Gf, b, m)
    b = reshape(b, m, m, m);
    vft = fftn(b,size(Gf));
    vft = vft.*Gf;
    vft = ifftn(vft);
    v = vft(1:m,1:m,1:m);
    v = v(:);
end
