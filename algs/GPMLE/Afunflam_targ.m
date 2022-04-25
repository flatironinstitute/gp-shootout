function A = Afunflam_targ(i,j,x,ker,nsrc)
    m = length(i);
    n = length(j);
    A = zeros(m,n);
    A(i>nsrc,j<=nsrc) = densekermat(ker.k,x(:,i(i>nsrc)),x(:,j(j<=nsrc)));
end