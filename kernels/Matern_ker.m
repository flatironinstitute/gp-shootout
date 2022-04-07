function ker = Matern_ker(dim, nu, l)
% MATERN12_KER   Matern nu=1/2 kernel func and Fourier transform, general dim
%
% ker = Matern12_ker(dim,l) returns ker.k a function handle taking an array of
%  distances to kernel values, and ker.khat, the function handle to its
%  Fourier transform (mapping wavevector magnitudes to values).
%  Latter uses Phillip's FT convention. l is the distance scale, aka rho
%
% To test see: TEST_KERS

if dim<1, error('dim must be at least 1!'); end
if l<=0, error('l lengthscale must be positive!'); end
if ~ismember(nu, [0.5, 1.5, 2.5]), error('for now, nu must be in [0.5, 1.5, 2.5]'); end 

if nu == 0.5
    ker.k = @(d) exp( (-1/l) * d);                                % d = |x|
elseif nu == 1.5
    ker.k = @(d) (1 + sqrt(3) .* d ./ l) .* exp(-sqrt(3) .* d ./ l);
elseif nu == 2.5
    ker.k = @(d) (1 + sqrt(5).*d./l + 5.*d.^2/(3*l^2)) .* exp(-sqrt(5) .* d / l);
end

% following taken from philip nufft_gps2d:  
scaling = (2*sqrt(pi))^dim * gamma(nu+dim/2) * (2*nu)^nu / (gamma(nu) * l^(2*nu));   % yuk; see https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
ker.khat = @(xid) scaling * (2*nu/l^2 + (4*pi^2) * xid.^2).^(-(nu + dim/2));
% xid means |xi|
