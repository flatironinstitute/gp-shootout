function [Kpxy,nbr] = pxyfunflam_targ(rc,xtrg,x,slf,nbr,l,ctr,proxy,shift,ker)
  pxy = proxy + shift.*l + ctr;  % scale and translate reference points
  if rc == 'r'
      Kpxy = densekermat(ker.k,xtrg(:,slf),pxy);
      dr = x(:,nbr)-ctr;
  else
      Kpxy = densekermat(ker.k,pxy,x(:,slf));
      dr = xtrg(:,nbr)-ctr;
  end
  [ndim,~] = size(dr);
  if(ndim == 1) 
    nbr = nbr(abs(dr)/l < 1.5);
  else
      
    nbr = nbr(sum((dr./l).^2) < 1.5^2);
  end
end
