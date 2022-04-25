function [Kpxy,nbr] = pxyfunflam_targ(x,slf,nbr,l,ctr,proxy,ker,nsrc)
  pxy = proxy.*l + ctr;  % scale and translate reference points
  [~,m] = size(pxy);
  n = length(slf);
  Kpxy = zeros(m,n);
  Kpxy(:,slf<=nsrc) = densekermat(ker.k,pxy,x(:,slf(slf<=nsrc)));
  % proxy points form ellipse of scaled "radius" 1.5 around current box
  % keep among neighbors only those within ellipse
  nbr = nbr(sum(((x(:,nbr) - ctr)./l).^2) < 1.5^2);
end