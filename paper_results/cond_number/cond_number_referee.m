dim = 1;
l = 0.1;
ker = Matern_ker(dim, 0.5, l);

rng(1)
ndouble = 4;

sig = 0.3;
sigmasq = sig*sig;

eps = 1e-3;

x_i = rand(100000,1);
y_i = rand(100000,1);
x_i(1) = 0;
x_i(2) = 1;

figure(1)
clf

figure(2)
clf

for i=1:ndouble
    n = 100*10^(i-1);

    [beta, xis, ~, ~, A, X, ws] = efgp1d_dense(x_i(1:n), y_i(1:n), ...
        sigmasq, ker, eps, x_i(1:n));
    
    [n1,~] = size(A);
    % A = A + sigmasq*eye(n1);
    [s] = svd(A);
    
    figure(1)
    semilogy((1:n1), s/s(1), '.'); hold on;

    s2 = s + sigmasq;
    figure(2)
    semilogy((1:n1), s2/s2(1), '.'); hold on;
end

%%
ntol = 8;
nn_use = 2;
niters_cg = zeros(ntol,nn_use);
niters_gmres = zeros(ntol,nn_use);


for j=1:nn_use
    n = 10000*10^(j-1);

    [beta, xis, ~, ~, A, X, ws] = efgp1d_dense(x_i(1:n), y_i(1:n), ...
            sigmasq, ker, eps, x_i(1:n));

    A = A + sigmasq*eye(n1);

    rhs = X'*y_i(1:n);
    maxit = n1*10;
    for i=1:ntol
        cgtol = 10^(-i);
        [beta,flag,relres,nn,resvec] = gmres(A, rhs, [], cgtol, maxit);  
        niters_gmres(i,j) = nn(2);
        [beta,flag,relres,niters_cg(i,j),resvec] = pcg(A, rhs, cgtol, maxit);  
    end
end

figure(3)
semilogx(10.^(-(1:ntol)).', niters_cg(:,1), 'k.', 'MarkerSize', 10); hold on;
semilogx(10.^(-(1:ntol)).', niters_gmres(:,1), 'r.', 'MarkerSize', 10);
semilogx(10.^(-(1:ntol)).', niters_cg(:,2), 'ks', 'MarkerSize', 10); 
semilogx(10.^(-(1:ntol)).', niters_gmres(:,2), 'rs', 'MarkerSize', 10);
