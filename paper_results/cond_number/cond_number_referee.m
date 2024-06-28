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

for i=1:ndouble
    n = 100*10^(i-1);

    [beta, xis, ~, ~, A, X, ws] = efgp1d_dense(x_i(1:n), y_i(1:n), ...
        sigmasq, ker, eps, x_i(1:n));
    
    [n1,~] = size(A);
    % A = A + sigmasq*eye(n1);
    [s] = svd(A);
    
    s2 = s + sigmasq;
    figure(1)
    legendi = sprintf('N = %0.0e', n);
    semilogy((1:n1), s2/s2(1), '.', 'DisplayName', legendi); hold on;
end

figure(1)
legend
title('normalized eigenvalues \lambda_j of X^*X + \sigma^2 for various N', 'FontSize', 16)
xlabel('j', 'FontSize', 16) 
ylabel('\lambda_j', 'FontSize', 16) 


filename = 'spectra.pdf';
exportgraphics(figure(1), filename)