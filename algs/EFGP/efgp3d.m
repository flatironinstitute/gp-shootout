function [beta, xis, ytrg, iter, time_info] = efgp3d(x, y, sigmasq, ker, eps, xsol,opts)
% EFGP3D   equispaced Fourier NUFFT-based GP regression in 3D
%
% [beta, xis, yhat, iter, time_info] = efgp3d(x, y, sigmasq, ker, eps, xsol)
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
% opts - [optional] struct controlling method params including:
%         cg_tol_fac - cg tolerance will be set to tol/cg_tol_fac, default
%         value 1;
%         use_integral - whether to use old tail integral estimates for
%         determining h,m
%         l2scaled - whether to use l2 scaling of integral with new
%         heuristics in determining h,m
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

  if(nargin == 6), opts = []; end
  tic_precomp = tic;
  x0 = min(x); x1 = max(x);   % both row 2-vectors
  L = max(x1-x0);    % worst-axis domain length *** could check xsol too?
  
  [xis, h, mtot] = get_xis(ker, eps, L,opts);   % mtot = 2m+1 from paper
  [xis_xx, xis_yy, xis_zz] = ndgrid(xis,xis,xis);        % assumes isotropic
  % center all coords for NUFFTs domain, then do 2pi.h ("tph") rescaling...
  xcen = (x1+x0)/2;                    % row vec
  tphx = 2*pi*h*(x - xcen);            % note broadcast over rows
  tphxsol = 2*pi*h*(xsol - xcen);      % "

  % khat & quadr weights of Fourier basis funcs
  rs = sqrt(xis_xx.^2 + xis_yy.^2 + xis_zz.^2);
  dim = 3; ws = sqrt(ker.khat(rs) * h^dim);
    
  % precomputation for fast apply of X*X
  nuffttol = eps/10;
  Gf = getGf3d(nuffttol, tphx, mtot);
  
  % conjugate gradient
  ws_flat = ws(:);
  Afun = @(a) ws_flat .* apply_xtx3d(Gf, ws_flat .* a, mtot) + sigmasq .* a;
    
  isign = -1;
  rhs = finufft3d1(tphx(:,1),tphx(:,2),tphx(:,3), y, isign, nuffttol, mtot, mtot, mtot);
  rhs = reshape(rhs .* ws, [], 1);
  t_precomp = toc(tic_precomp);
    
    % solve linear system
    tic_cg = tic; 
    cg_tol_fac = 1;
    if(isfield(opts,'cg_tol_fac'))
        cg_tol_fac = opts.cg_tol_fac;
    end
    cgtol = eps/cg_tol_fac;
    maxit = 3*mtot^3;
    [beta,flag,relres,iter,resvec] = pcg(Afun, rhs, cgtol, maxit);
    t_cg = toc(tic_cg);

    % evaluate posterior mean 
    tic_post = tic;
    wsbeta = ws .* reshape(beta, [mtot, mtot, mtot]);
    isign = +1;
    yhat = finufft3d2(tphxsol(:,1),tphxsol(:,2),tphxsol(:,3), isign, nuffttol, wsbeta);
    ytrg.mean = real(yhat);
    t_post = toc(tic_post);

    % package times for output
    time_info = [t_precomp, t_cg, t_post];
    time_info = [time_info, sum(time_info)];
end


function [Gf] = getGf3d(nuffttol, tphx, mtot)
% array to multiply by on the 3D FFT side that convolves with Toeplitz col blk
    N = size(tphx,1);
    c = complex(ones(N, 1));            % unit strengths
    isign = -1;
    o.modeord = 1;
    XtXcol_blk = finufft3d1(tphx(:,1),tphx(:,2),tphx(:,3), c, isign, nuffttol, 2*mtot-1,2*mtot-1,2*mtot-1, o);
    Gf = fftn(XtXcol_blk);
end


function [v] = apply_xtx3d(Gf, b, mtot)
    b = reshape(b, mtot, mtot, mtot);
    vft = fftn(b,size(Gf));
    vft = vft.*Gf;
    vft = ifftn(vft);
    v = vft(1:mtot,1:mtot,1:mtot);
    v = v(:);
end
