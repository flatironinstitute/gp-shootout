function [Kpxy,nbr] = pxyfunflam(x,slf,nbr,l,ctr,proxy,shift,ker)
  pxy = proxy +shift.*l + ctr;  % scale and translate reference points
  Kpxy = densekermat(ker.k,pxy,x(:,slf));
  % proxy points form ellipse of scaled "radius" 1.5 around current box
  % keep among neighbors only those within ellipse
  [ndim,~] = size(pxy);
  if(ndim==1)
      nbr = nbr(abs(x(:,nbr) - ctr)./l < 1.5);
  else
    nbr = nbr(sum(((x(:,nbr) - ctr)./l).^2) < 1.5^2);
  end
end
