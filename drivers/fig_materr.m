% script to make figures showing good heuristic Matern truncation error estimate
% Splits out computation into a func.
% Barnett 8/31/22 based on Rachh.
clear

eps = 1e-6;  % targ aliasing err to control h
var = 1;  % always

figure;
for dim=[1 2 3]
  fprintf('dim=%d...\n',dim); subplot(1,3,dim);
  nuu = [0.5 2.5]; %[0.5 1.5 2.5];  % nu values to test
  colors = 'kr';
  for j=1:numel(nuu), nu = nuu(j); color = colors(j);
    ll = [0.02 0.2];   % ell values to test
    marks = 'o+';
    for i=1:numel(ll), l = ll(i); mark = marks(i);
      
      ker = Matern_ker(dim,nu,l,var);

      h = 1/(1+0.7*l/sqrt(nu)*log(1/eps));   % alex ok heuristic for aliasing
  
      aliaserr = ker.k(1/h-1);   % nearest image (good unif err estim)
      fprintf('\tnu=%.3g\tl=%.3g\th=%.3g\t(est max alias err=%.3g)\n',nu,l,h,aliaserr)
      
      % go up to N=1e6 modes total
      if dim==1, ifac = 1:0.5:6;
      elseif dim==2, ifac = 1:0.25:3;
      elseif dim==3, ifac = 0.6:0.2:2;
      end
      mms = floor(10.^(ifac));
      
      errs = materr_vs_m(ker,mms, h, 1e-2*eps);    % meas L2 error (NUFFT+quadr)
      ph(i,j) = loglog(mms,errs,[color mark '-'],'MarkerSize',10); hold on;
      legstr{i,j} = sprintf('$\\nu=%g, \\ell=%g$',nu,l);
      
      rfac = nu^(nu-1)/2^(nu)/pi^(dim/2+2*nu)*gamma(nu+0.5)/gamma(nu);
      errs_heur = rfac ./ (h*mms).^(2*nu+dim/2) ./ l^(2*nu);
      loglog(mms,errs_heur,[color mark '--'],'MarkerSize',5);
      %title(sprintf('test Matern ker accuracy: d=%d, nu=%g, ell=%g, h=%g',dim,nu,l,h))
      %hline(aliaserr,'g','estim alias err')
    end
  end
  xlabel('m');
  ylabel('RMS error in $\tilde{k}$','interpreter','latex');
  axis tight; v = axis; axis([v(1:2) 0.1*eps 1]);   % clip vertical domain
  legend(ph(:),legstr{:},'interpreter','latex');
  title(sprintf('Matern, d=%d',dim))
end

%%%%%%%%%%%%
function errs = materr_vs_m(ker,mms, h, epsnufft)
  % use quadrature to estimate expected rms Matern matrix element approx err
  % via pointy-hat weighted single integral over [-1,1]^d, which in turn equals
  % the double integral of |kapprox - k|^2 over [0,1]^d cross itself.
  %  Inputs: ker = kernel struct
  %          mms = list of m values to test
  %          h = quad spacing in xi, using exp(2pi.i.x.xi) convention.
  %          epsnufft = NUFFT tolerance

  dim = ker.dim; nu =ker.nu; l = ker.l; % get stuff out of ker object
  if dim==1      % we guess enough to capture small-l behavior
    nleg = 500;
  elseif dim==2
    nleg = 100;
  elseif dim==3
    nleg = 30;
  end
  % seems to build spectral 2-panel quadr over [-1,1]
  [xleg,wleg] = legpts(nleg);
  xleguse = [-1+(xleg+1)/2; (xleg+1)/2];
  wleguse = [wleg/2, wleg/2];           % row vec
  wleguse = (1 - abs(xleguse.')).*wleguse;  % top-hat weight: conv(rho,rho)
  
  if(dim == 1)
    xq = xleguse(:);
    wq = wleguse(:);
    xuse = xq*2*pi*h;
    rr = abs(xq);
  elseif dim==2
    [x,y] = meshgrid(xleguse);
    [wy,wx] = meshgrid(wleguse);
    xq = x(:);
    yq = y(:);
    wq = wx(:).*wy(:);
    xuse = xq*2*pi*h;        % rescaled space targets in exp(ikx) convention
    yuse = yq*2*pi*h;
    rr = sqrt(xq.^2 + yq.^2);    
  elseif dim==3
    [x,y,z] = ndgrid(xleguse,xleguse,xleguse);
    [wx,wy,wz] = ndgrid(wleguse,wleguse,wleguse);
    xq = x(:); yq = y(:); zq = z(:);
    wq = wx(:).*wy(:).*wz(:);
    xuse = xq*2*pi*h;        % rescaled space targets in exp(ikx) convention
    yuse = yq*2*pi*h;
    zuse = zq*2*pi*h;
    rr = sqrt(xq.^2 + yq.^2 + zq.^2);  
  end
  fprintf('\ttot real-space quad nodes = %d\n',numel(xuse))

  
  %fprintf('pointy-hat wei quadr over [-1,1]^d test error: %.3g\n',abs(sum(wq)-1))
  
  kvals = ker.k(rr);

  errs = zeros(size(mms));
  for i=1:length(mms)
    mmax = mms(i);
    xis = -mmax:1:mmax;
    xis = xis*h;           % Fourier regular nodes
    
    if(dim == 1)
      xisr = abs(xis);
    elseif(dim == 2)
      [xisy,xisx] = meshgrid(xis);
        xisr = sqrt(xisx.^2 + xisy.^2);
    elseif(dim == 3)
      [xisx,xisy,xisz] = ndgrid(xis,xis,xis);
      xisr = sqrt(xisx.^2 + xisy.^2 + xisz.^2);
    end

    % Now evaluate truncated fourier approximation
    fvals = ker.khat(xisr);

    epsfudge = 1e-2;
    if(dim == 1)
      cvals = h^(dim)*finufft1d2(xuse,1,epsnufft,fvals);
    elseif(dim == 2)
      cvals = h^(dim)*finufft2d2(xuse,yuse,1,epsnufft,fvals);
    elseif(dim == 3)
      cvals = h^(dim)*finufft3d2(xuse,yuse,zuse,1,epsnufft,fvals);
    end

    errs(i) = sqrt(sum((kvals-real(cvals)).^2.*wq));
end
end  %%%%%%%%%%

