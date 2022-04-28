function A = Afunflam_targ(i,j,xtrg,x,ker)
    A = densekermat(ker.k,xtrg(:,i),x(:,j));
end