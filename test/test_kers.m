% script to test all kernel functions have roughly the right Fourier transform,
% using equispaced quadrature (PTR) Kronecker product over dims.
% For now, only SE kernel and Matern-1/2, Matern-3/2, Matern-5/2.
% *** to do: more formal unit tests
% Barnett 4/5/22 started.
clear; disp('test_kers: user should visually check that the below errors are not O(1). This checks that prefactors are not O(1) wrong...');

l = 0.1;     % spatial scale (for all kernels that have one)
xid = 3.0;   % freq mag |xi| to test (let's make all kernels use same for now)

for type=1:4   % -------------- kernel types --------------
  for dim=1:3      % .......
    switch type        
      case 1
        nam = 'SE';
        ker = SE_ker(dim,l);
        h = l/2;         % small enough spacing for emach on Gaussian
        L = l*sqrt(2*log(1/eps));   % max x to get tails to emach
      case 2
        nam = 'Matern nu=1/2';
        nu = 0.5;
        ker = Matern_ker(dim, nu, l);
        h = l/5;    % crappy convergence
        tol = 1e-5; L = l*log(1/tol);   % max x to get tails to err tol
      case 3
        nam = 'Matern nu=3/2';
        nu = 1.5;
        ker = Matern_ker(dim, nu, l);
        h = l/5;    % less crappy convergence
        tol = 1e-5; L = l*log(1/tol);   % max x to get tails to err tol
      case 4
        nam = 'Matern nu=5/2';
        nu = 2.5;
        ker = Matern_ker(dim, nu, l);
        h = l/5;    % even less crappy convergence
        tol = 1e-5; L = l*log(1/tol);   % max x to get tails to err tol
      otherwise
        error(sprintf('ker type %d not implemented in test_kers!',type));
    end
    khat = ker.khat(xid);      % get analytic FT claim
    n = ceil(2*L/h);
    x = -L+(0:n-1)*h;
    if dim==1
      fprintf('type=%d: %s...\n',type,nam)
      khatap = h*sum(ker.k(abs(x)).*cos(2*pi*xid*x));
    elseif dim==2 
      [xx yy] = ndgrid(x,x);
      n = numel(xx);
      dd = sqrt(xx(:).^2+yy(:).^2);        % list of distances
      khatap = h^2*sum(ker.k(dd).*cos(2*pi*xid*xx(:)));   % xi = (xid,0)
    elseif dim==3
      [xx yy zz] = ndgrid(x,x,x);
      n = numel(xx);
      dd = sqrt(xx(:).^2+yy(:).^2+zz(:).^2);        % list of distances
      khatap = h^3*sum(ker.k(dd).*cos(2*pi*xid*xx(:)));   % xi = (xid,0,0)
    end
    fprintf('dim=%d (%d quadr pts):    \tkhat rel err = %.3g\n',dim,n,(khatap-khat)/khat)
  end         % .........
end         % ------------------
