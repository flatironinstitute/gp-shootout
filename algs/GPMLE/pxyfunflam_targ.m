function [Kpxy,nbr] = pxyfunflam_targ(rc,xtrg,x,slf,nbr,l,ctr,proxy,ker)
  pxy = proxy.*l + ctr;  % scale and translate reference points
  if rc == 'r'
      Kpxy = densekermat(ker.k,xtrg(:,slf),pxy);
      dr = x(:,nbr)-ctr;
  else
      Kpxy = densekermat(ker.k,pxy,x(:,slf));
      dr = xtrg(:,nbr)-ctr;
  end
  nbr = nbr(sum((dr./l).^2) < 1.5^2);
end