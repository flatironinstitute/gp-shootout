dim = 1;
l = 0.1;
ker = Matern_ker(dim, 0.5, l);

rng(1)
ndouble = 5;

for i=1:ndouble
    n = 200*2^(i);
    x_i = rand(2,n);
    K = densekermat(ker.k, x_i);

    [s] = svd(K);
    
    semilogy((1:n)/n, s/s(1), '.'); hold on;
end
