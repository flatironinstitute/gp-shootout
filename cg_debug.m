% data
N = 1e5;
f = @(x) cos(6*2*pi*x) / 2;
sigmatrue = 0.5;
dim = 1;
[x, meas, truemeas] = get_randdata(dim, N, f, sigmatrue);

% targets
ntrgs = 100;
xtrgs = linspace(0, 1, ntrgs);

% kernel
sigmasq = sigmatrue^2;
l = 0.1;
dim = 1;
ker = SE_ker(dim, l);

opts.tol = 1e-10;
[ytrgs, info] = EFGP_old(x, meas, sigmasq, ker, xtrgs, opts);
[ytrgs, info] = EFGP_old_dense(x, meas, sigmasq, ker, xtrgs, opts);
opts.tol = 1e-12;
[ytrgs2, info2] = EFGP_old(x, meas, sigmasq, ker, xtrgs, opts);
[ytrgs2, info2] = EFGP_old_dense(x, meas, sigmasq, ker, xtrgs, opts);

rms_err = rms(ytrgs.mean - ytrgs2.mean);
m = size(info.xis);
fprintf('m: %g \n', m);
fprintf('iters: %g \n', info.iter);
fprintf('rms err: %g \n', rms_err);




%%%%%%%%%%%%%%%%%%%%%%%
%%%% old functions %%%%
%%%%%%%%%%%%%%%%%%%%%%%


function [ytrg, info] = EFGP_old(x, meas, sigmasq, ker, xtrg, opts)
    xsol = xtrg';          % and transpose to Philip n*d shape   

    [info.beta, info.xis, yhat, info.iter, cpu_time] = function_space1d(x', meas, sigmasq, ker.khat, ker.k, opts.tol, xsol);

    info.h = info.xis(2)-info.xis(1); info.ximax = max(info.xis);
    info.cpu_time.precomp = cpu_time(1);
    info.cpu_time.cg = cpu_time(2);
    info.cpu_time.mean = cpu_time(3);
    ytrg.mean = yhat;
end


function [beta, xis, yhat, iter, ts] = function_space1d(x, y, sigma2, khat, k, eps, xsol)
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
    tol = eps / 100;
    c = ones(N, 1);
    XtXrow = finufft1d1(x*2*pi*h,c,+1,tol,2*m-1)'; 
    Gf = fftn(XtXrow.');

    % construct rhs with fft
    isign = -1;
    tol = eps / 100;
    rhs = finufft1d1(2*pi*x*h,y,isign,tol,m);
    rhs = ws .* rhs;
    
    t_precomp = toc(tic_precomp);
    
    % solve linear system with conjugate gradient
    Afun = @(a) ws .* Afun2(Gf, ws .* a) + sigma2 .* a; 
    
    tic_cg = tic; 
    itol = eps; 
    [beta,flag,relres,iter,resvec] = pcg(Afun,rhs,itol, m);
    t_cg = toc(tic_cg);
    
    % tabulate solution using fft
    tol = eps / 100;
    tmpvec = ws .* beta;
    tic_post = tic;
    yhat = finufft1d2(2*pi*h*xsol,+1,tol, tmpvec);
    yhat = real(yhat);
    t_post = toc(tic_post);

    ts = [t_precomp, t_cg, t_post];

end


function [v2] = Afun2(Gf, a)
    m = numel(a);
    af = fftn(a, size(Gf));
    vft = af .* Gf;
    vft = ifftn(vft);
    v2 = vft(m:end);
end


function [ytrg, info] = EFGP_old_dense(x, meas, sigmasq, ker, xtrg, opts)
    xsol = xtrg';         
    [info.beta, info.xis, yhat, info.A, iter] = function_space1d_naive(x', meas, sigmasq, ker.khat, ker.k, opts.tol, xsol);
    info.h = info.xis(2)-info.xis(1); info.ximax = max(info.xis);
    info.iter = iter;
    ytrg.mean = yhat;
end


function [beta, xis, yhat, A, iter] = function_space1d_naive(x, y, sigma2, khat, k, eps, xsol)
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
    ws = sqrt(khat(xis) * h);
    X = exp(1i * x * xis * 2 * pi);
    A = diag(ws) * (X'* X) * diag(ws) + sigma2 * eye(m);
    
    % solve A beta = 
    beta = A \ (diag(ws) * X' * y);





    % solve linear system with conjugate gradient
    Afun = @(a) ws .* Afun2(Gf, ws .* a) + sigma2 .* a; 
    Afun = @(a) A * a;
    rhs = diag(ws) * X' * y;
    tic_cg = tic; 
    itol = eps; 
    [beta,flag,relres,iter,resvec] = pcg(Afun,rhs,itol, 2*m);
    t_cg = toc(tic_cg);





    % tabulate posterior mean
    Xsol = exp(1i * xsol * xis * 2 * pi) * diag(ws);
    yhat = Xsol * beta;
    yhat = real(yhat);
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