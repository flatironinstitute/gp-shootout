function [Kpxy,nbr] = pxyfunflam(x,slf,nbr,l,ctr,proxy,ker)
  pxy = proxy.*l + ctr;  % scale and translate reference points
  Kpxy = densekermat(ker.k,pxy,x(:,slf));
  % proxy points form ellipse of scaled "radius" 1.5 around current box
  % keep among neighbors only those within ellipse
  nbr = nbr(sum(((x(:,nbr) - ctr)./l).^2) < 1.5^2);
end