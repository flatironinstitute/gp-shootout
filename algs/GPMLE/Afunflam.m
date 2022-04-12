function A = Afunflam(i,j,x,ker,sigmasq)
    A = densekermat(ker.k,x(:,i),x(:,j));
    [I,J] = ndgrid(i,j);
    idx = I == J;
    A(idx) = A(idx) + sigmasq;
end