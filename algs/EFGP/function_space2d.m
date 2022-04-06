function [beta, xis, yhat, iter, time_info] = function_space2d(eps, ker, xs, y, sigma2, xsols)
%  *** to rationalize interface
%  *** to doc
%  *** to unify and split out getL below
%
  
% get kernel functions
  k = ker.k; khat = ker.khat;

    % support of functionin time domain
    tic_precomp = tic;
    Ltime = getL(eps, k);
    Ltime = max(1, Ltime);
    
    % nyquist
    hnyq = 1/(2*Ltime);
    
    % m must be odd!
    Lfreq = getL(eps, khat);

    % discretization in Fourier domain
    m = 2*Lfreq/hnyq + 1;
    m = 2*ceil(m/2) + 1;
    xis = linspace(-Lfreq, Lfreq, m);
    h = xis(2) - xis(1);
    %%%fprintf('number of basis functions =\n %.3g\n', m^2)
    [xis_xx, xis_yy] = meshgrid(xis);

    const = h;
    rs = sqrt(xis_xx.^2 + xis_yy.^2);
    ws = sqrt(khat(rs)) * const;
    
    % precomputation for fast apply of X*X
    Gf = getGf(xs, xis);
    
    % conjugate gradient
    ws_flat = reshape(ws, m^2, 1);
    Afun = @(a) ws_flat .* apply_xtx(Gf, ws_flat .* a) + sigma2 .* a;
    
    isign = -1;
    eps1 = eps / 100;
    rhs = finufft2d1(2*pi*xs(:,1)*h, 2*pi*xs(:,2)*h, y, isign, eps1, m, m);
    rhs = reshape(rhs.' .* ws, [], 1); 
    
    t_precomp = toc(tic_precomp);
    % solve linear system
    tic_cg = tic;
    [beta,flag,relres,iter,resvec] = pcg(Afun,rhs,eps1,m^2);
    t_cg = toc(tic_cg);
    
    % evaluate posterior mean 
    tic_post = tic;
    tmpvec = ws .* reshape(beta, [m, m]);
    tmpvec = tmpvec.';
    isign = +1;
    yhat = finufft2d2(2*pi*h * xsols(:,1), 2*pi*h * xsols(:,2),isign,eps1,tmpvec);
    t_post = toc(tic_post);

    % package timings
    time_info = [t_precomp, t_cg, t_post];

    % convert to real
    yhat = real(yhat);

end


function L = getL(eps, kern_hat)
% find numerical bandlimit of spectral density using bisection
    k0 = kern_hat(0);
    a = 0;
    b = 1000;
    nmax = 10;
    for i=1:nmax
        if (kern_hat(b)/k0) > eps
            b = 2*b;
        else
            break
        end
    end
    %%%fprintf('starting upper bound increased to: \n %.3g\n', b);

    % start bisection
    nmax = 200;
    for i=1:nmax
        mid = (a + b)/2;
        fmid = kern_hat(mid)/k0;
        if fmid > eps
            a = mid;
        else
            b = mid;
        end
    end
    L = mid;
end


function [Gf] = getGf(xs, z)
    N = length(xs);
    m = length(z);
    % determine all the differences where kernel is to be evaluated...
    % ...for both x and y directions 
    h = z(2) - z(1);

    % parameters for fft
    c = 0i + ones(N, 1);
    isign = -1;
    eps = 1e-10;
    
    o.modeord = 1;
    Gf = finufft2d1(2*pi*xs(:,1)*h,2*pi*xs(:,2)*h,c,isign,eps,2*m-1, 2*m-1, o);
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

