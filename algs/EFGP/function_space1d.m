function [beta, xis, yhat, iter, time_info] = function_space1d(x, y, sigma2, ker, eps, xsol)
% *** to doc
%
% Outputs:
% beta - vector of Fourier basis weights (not really for the user)
% xis  - Fourier freqs used (not really for the user)
% yhat - posterior means at xsol ordinates   <- the only user output
% iter - diagnostics from CG
% ts   - diagnostic list of timings
    
  k = ker.k; khat = ker.khat;    % get functions, new kernel format
  N = numel(y);

    tic_precomp = tic;
    % find support of function in time domain
    Ltime = getL(eps, k);
    Ltime = max(1, Ltime);

    % nyquist
    hnyq = 1/(2*Ltime);
    
    % find numerical bandlimit via bisection
    Lfreq = getL(eps, khat);
    
    % number of nodes for nyquist, should be odd
    m = 2*Lfreq/hnyq + 1;
    m = 2*ceil(m/2) + 1;
    xis = linspace(-Lfreq, Lfreq, m);
    h = xis(2) - xis(1);
    
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
    Afun = @(a) ws .* Afun2(Gf, ws .* a) + sigma2 .* a; 
    
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
end


function L = getL(eps, kern_hat)
% find numerical bandlimit of spectral density using bisection
    a = 0;
    b = 1000;
    nmax = 10;
    for i=1:nmax
        if kern_hat(b) > eps
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
        fmid = kern_hat(mid);
        if fmid > eps
            a = mid;
        else
            b = mid;
        end
    end
    L = mid;
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