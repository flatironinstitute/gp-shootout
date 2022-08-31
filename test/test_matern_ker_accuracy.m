% script to estimate L2([-1,1]^d) norm of error tilde{k}(x)-k(x) for EFGP,
% compare to heuristic over range of m.
% Manas Rachh; Alex Barnett attempting nu.neq.0.5 improvements 8/31/22

clear
nu = 2.5;
l = 0.01;   % try down to 0.03, smaller = better for trunc err assessment
dim = 1;
var = 1;  % always
ker = Matern_ker(dim,nu,l,var);
eps = 1e-8;
%h = 1/(1+l*sqrt(2*dim/nu)*log(dim*3^dim/eps));   % fixed from Cor 9.
h = 1/(1+0.7*l/sqrt(nu)*log(1/eps));   % alex ok heuristic for aliasing
%h = h/2;            % why? make sure h-converged, so trunc error dominates.

aliaserr = ker.k(1/h-1);   % nearest image (good unif err estim)
fprintf('h=%.3g, estim max aliasing err = %.3g\n',h,aliaserr)

nleg = 300;         % seems to build spectral 2-panel quadr over [-1,1]
[xleg,wleg] = legpts(nleg);
xleguse = [-1+(xleg+1)/2; (xleg+1)/2];
wleguse = [wleg/2, wleg/2];
wleguse = (1 - abs(xleguse.')).*wleguse;  % apply top-hat weight: conv(rho,rho)

if(dim == 2)

    

    [x,y] = meshgrid(xleguse);
    [wy,wx] = meshgrid(wleguse);

    xq = x(:);
    yq = y(:);
    wq = wx(:).*wy(:);
    xuse = xq*2*pi*h;        % rescaled space targets in exp(ikx) convention
    yuse = yq*2*pi*h;
    rr = sqrt(xq.^2 + yq.^2);

elseif(dim == 1)
   xq = xleguse(:);
   wq = wleguse(:);
   xuse = xq*2*pi*h;
   rr = abs(xq);
end

fprintf('tri-hat wei quadr over [-1,1]^d test error: %.3g\n',abs(sum(wq)-1))

kvals = ker.k(rr);

if dim==1
  ifac = 1:0.25:5;
elseif dim==2
  ifac = 1:0.25:3;
elseif dim==3
  ifac = 0.6:0.2:2;
end
mms = floor(10.^(ifac));

errs = zeros(size(mms));
for i=1:length(mms)
    mmax = mms(i);
    xis = -mmax:1:mmax;
    xis = xis*h;           % Fourier regular nodes
    
    if(dim == 2)
        [xisy,xisx] = meshgrid(xis);
        xisr = sqrt(xisx.^2 + xisy.^2);
    elseif(dim == 1)
        xisr = abs(xis);
    end
    




    % Now evaluate truncated fourier approximation
    fvals = ker.khat(xisr);
    if(dim == 2)
        cvals = h^(dim)*finufft2d2(xuse,yuse,1,1e-13,fvals);
    else
        cvals = h^(dim)*finufft1d2(xuse,1,1e-13,fvals);
    end


    errs(i) = sqrt(sum((kvals-real(cvals)).^2.*wq));
end

% 
% kvals = reshape(kvals,size(x));
% cvals2 = reshape(cvals2,size(x));

% figure(1)
% clf()
% hh = pcolor(x,y,kvals); 
% set(hh,'EdgeColor','none');    
% colorbar();
% 
% figure(2)
% clf()
% hh = pcolor(x,y,abs(cvals2)); 
% set(hh,'EdgeColor','none');    
% colorb/ar();

figure(3)
clf
loglog(mms,errs,'k.-','MarkerSize',20); hold on;
rfac = nu^(nu-1)/2^(nu)/pi^(dim/2+2*nu)*gamma(nu+0.5)/gamma(nu);
%errs_ex = rfac./mms.^(2*nu+dim/4)/h.^(2*nu+(dim-1)/2)/l;  % too small for nu>=3/2
errs_ex = rfac ./ (h*mms).^(2*nu+dim/4) ./ l^(2*nu);    % skip h for now?
loglog(mms,errs_ex,'r.-','MarkerSize',20)
legend('estimated error','heuristic');
title(sprintf('test Matern ker accuracy: d=%d, nu=%g, ell=%g, h=%g',dim,nu,l,h))
xlabel('m'); ylabel('$L^2([-1,1]^d)$ error of $\tilde{k}$','interpreter','latex');
hline(aliaserr,'g','estim alias err')
