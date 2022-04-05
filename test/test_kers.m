% script to test all kernel functions have roughly the right Fourier transform,
% using equispaced quadrature (PTR) Kronecker product over dims.
% For now, only SE kernel.
% Barnett 4/5/22
clear

xid = 3.0; % freq mag |xi| to test (let's make all kernels use same for now)
for type=1:1   % -------------- kernel types --------------
  for dim=1:3      % .......
    switch type
      case 1        % SE ker
        l = 0.1;       % spatial scale
        ker = SE_ker(dim,l);
        h = l/2;         % small enough spacing for emach on Gaussian
        L = l*sqrt(2*log(1/eps));   % max x to get tails to emach
      otherwise
        error('ker type not implemented in test!');
    end
    khat = ker.khat(xid);      % get analytic FT claim
    x = -L+(0:2*L/h)*h;
    if dim==1
      khatap = h*sum(ker.k(abs(x)).*cos(2*pi*xid*x));
    elseif dim==2 
      [xx yy] = ndgrid(x,x);
      dd = sqrt(xx(:).^2+yy(:).^2);        % list of distances
      khatap = h^2*sum(ker.k(dd).*cos(2*pi*xid*xx(:)));   % xi = (xid,0)
    elseif dim==3
      [xx yy zz] = ndgrid(x,x,x);
      dd = sqrt(xx(:).^2+yy(:).^2+zz(:).^2);        % list of distances
      khatap = h^3*sum(ker.k(dd).*cos(2*pi*xid*xx(:)));   % xi = (xid,0,0)
    end
    fprintf('dim=%d:  khat rel err = %.3g\n',dim,(khatap-khat)/khat)
  end         % .........
end         % ------------------
