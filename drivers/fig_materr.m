% script to make figures showing good heuristic Matern truncation error estimate
% Splits out computation into a func.
% Barnett 8/31/22 based on Rachh.
clear

eps = 1e-8;  % targ aliasing err to control h
var = 1;  % always

errs_all = cell(3,3,2);
errs_rel_all = cell(3,3,2);

errs_heur_all = cell(3,3,2);
errs_rel_heur_all = cell(3,3,2);
rl2sq_all = cell(3,3,2);

figure(1);
clf
hold on;
figure(2);
clf
hold on;
for dim=[1 2 3]
  figure(1)  
  fprintf('dim=%d...\n',dim); subplot(1,3,dim);
  figure(2)  
  fprintf('dim=%d...\n',dim); subplot(1,3,dim);
  nuu = [0.5 1.5 2.5];  % nu values to test
  colors = 'krg';
  for j=1:numel(nuu), nu = nuu(j); color = colors(j);
    ll = [0.02 0.2];   % ell values to test
    marks = 'o+';
    for i=1:numel(ll), l = ll(i); mark = marks(i);
      
      ker = Matern_ker(dim,nu,l,var);

      h = 1/(1+0.85*l/sqrt(nu)*log(1/eps));   % alex ok heuristic for aliasing
  
      aliaserr = ker.k(1/h-1);   % nearest image (good unif err estim)
      fprintf('\tnu=%.3g\tl=%.3g\th=%.3g\t(est max alias err=%.3g)\n',nu,l,h,aliaserr)
      
      % go up to N=1e6 modes total
      if dim==1, ifac = 1:0.5:6;
      elseif dim==2, ifac = 1:0.25:3;  % for faster, only go to 3
      elseif dim==3, ifac = 0.8:0.2:1.8;   % for faster, only go to 1.8
      end
      mms = floor(10.^(ifac));
      
      [errs,errs_rel] = materr_vs_m(ker,mms, h, 1e-2*eps);    % meas L2 error (NUFFT+quadr)
      errs_all{dim,j,i} = errs;
      errs_rel_all{dim,j,i} = errs_rel;
      figure(1)
      ph(i,j) = loglog(mms,errs,[color mark '-'],'MarkerSize',10); hold on;
      legstr{i,j} = sprintf('$\\nu=%g,\\; \\ell=%g$',nu,l);  % note \ esc char

      % m-power:ok. l-power:ok. nu=1/2 heur needs prefac/=3.
      % nu=5/2 needs prefac*=3^{(d-1)/2} etc.
      %rfac = nu^(nu-1)/2^(nu)/pi^(dim/2+2*nu)*gamma(nu+0.5)/gamma(nu);
      rfac = 0.15/pi^(nu+dim/2);   %-(nu-1/2)*(dim-1)/4);  % no dim-dep
      errs_heur = rfac ./ (h*mms).^(2*nu+dim/2) ./ l^(2*nu);
      errs_heur_all{dim,j,i} = errs_heur;
      loglog(mms,errs_heur,[color mark '--'],'MarkerSize',5);
      
      figure(2)
      ph2(i,j) = loglog(mms,errs_rel,[color mark '-'],'MarkerSize',10); hold on;
      legstr2{i,j} = sprintf('$\\nu=%g,\\; \\ell=%g$',nu,l);  % note \ esc char
      
      rl2sq = (2*nu/pi/l^2)^(dim/2)*ker.khat(0)^2/2*gamma(dim/2+2*nu)/gamma(dim+2*nu)*2^(-dim/2);
      rl2sq_all{dim,j,i} = rl2sq;
      errs_heur = rfac/sqrt(rl2sq)./ (h*mms).^(2*nu+dim/2) ./ l^(2*nu);
      errs_rel_heur_all{dim,j,i} = errs_heur;
      loglog(mms,errs_heur,[color mark '--'],'MarkerSize',5);
      
      %title(sprintf('test Matern ker accuracy: d=%d, nu=%g, ell=%g, h=%g',dim,nu,l,h))
      %hline(aliaserr,'g','estim alias err')
    end
  end
  figure(1)
  xlabel('$m$','interpreter','latex');
  ylabel('RMS error in $\tilde{k}$','interpreter','latex');
  axis tight; v = axis; axis([v(1:2) 0.1*eps 1]);   % clip vertical domain
  legend(ph(:),legstr{:},'interpreter','latex');
  title(sprintf('Matern, d=%d',dim))
  
  figure(2)
  xlabel('$m$','interpreter','latex');
  ylabel('relative L2 error in $\tilde{k}$','interpreter','latex');
  axis tight; v = axis; axis([v(1:2) 0.1*eps 1]);   % clip vertical domain
  legend(ph2(:),legstr2{:},'interpreter','latex');
  title(sprintf('Matern, d=%d',dim))
end

set(gcf,'paperposition',[0 0 14 4])
print -dpng ../results/alex/fig_materr.png
% pdf is a pain


% check inversion formula to pick m:
m_heur = ( pi^(nu+dim/2)*l^(2*nu) * eps/0.15 )^(-1/(2*nu+dim/2)) / h;
err_heur = rfac ./ (h*m_heur).^(2*nu+dim/2) ./ l^(2*nu);
fprintf('check m_heuristic, should match: %g %g\n',eps,err_heur)



%%%%%%%%%%%%
function [errs_abs,errs_rel] = materr_vs_m(ker,mms, h, epsnufft)
  % use quadrature to estimate expected rms Matern matrix element approx err
  % via pointy-hat weighted single integral over [-1,1]^d, which in turn equals
  % the double integral of |kapprox - k|^2 over [0,1]^d cross itself.
  %  Inputs: ker = kernel struct
  %          mms = list of m values to test
  %          h = quad spacing in xi, using exp(2pi.i.x.xi) convention.
  %          epsnufft = NUFFT tolerance
  % Note: this repackages ../test/test_matern_ker_accuracy.m

  dim = ker.dim; nu =ker.nu; l = ker.l; % get stuff out of ker object
  if dim==1      % we guess enough to capture small-l behavior
    nleg = 500;
  elseif dim==2
    nleg = 100;
  elseif dim==3
    nleg = 50;
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

  errs_abs = zeros(size(mms));
  errs_rel = zeros(size(mms));
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

    
    if(dim == 1)
      cvals = h^(dim)*finufft1d2(xuse,1,epsnufft,fvals);
    elseif(dim == 2)
      cvals = h^(dim)*finufft2d2(xuse,yuse,1,epsnufft,fvals);
    elseif(dim == 3)
      cvals = h^(dim)*finufft3d2(xuse,yuse,zuse,1,epsnufft,fvals);
    end
    errs_abs(i) = sqrt(sum((kvals-real(cvals)).^2.*wq));
    errs_rel(i) = sqrt(sum((kvals-real(cvals)).^2.*wq))./sqrt(sum(kvals.^2.*wq));
end
end  %%%%%%%%%%

