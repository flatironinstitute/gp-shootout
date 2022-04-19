function ker = Matern_ker(dim, nu, l, var)
% MATERN_KER   Matern kernel func and Fourier transform, general nu and dim.
%
% ker = Matern_ker(dim,nu,l) returns the Matern kernel in general dimension dim,
%  with various nu parameters, and general lengthscale l (a.k.a. rho).
%  A struct is returned with at least the fields:
%    k -    a function handle taking an array of distances to kernel values.
%    khat - a function handle to its d-dim isotropic Fourier transform (mapping
%           wavevector magnitudes to values).
%    fam -  a descriptive string.
%    l, nu - lengthscale, smoothness param
%    dim  - spatial dimension.
%  nu can take one of only three values: 1/2, 3/2, or 5/2, currently.
%  See the paper for the Fourier transform convention (the 2pi is "upstairs").
%
% ker = Matern_ker(dim,nu,l,var) changes the variance (ie kernel value at r=0)
%  from its default of 1.0 (as above) to var.
%
% We follow the Matern definitions here:
%  https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
%  where they use sigma^2 for var, rho for l, and n for dim.
%
% To test see: TEST_KERS

if dim<1, error('dim must be at least 1!'); end
if l<=0, error('l lengthscale must be positive!'); end
if ~ismember(nu, [0.5, 1.5, 2.5]), error('for now, nu must be in [0.5, 1.5, 2.5]'); end 
if nargin<4, var = 1.0; end    % var is overall scale factor

if nu == 0.5
    ker.k = @(d) var*exp( (-1/l) * abs(d));                                % d = |x|
    ker.fam = 'matern12';
elseif nu == 1.5
    ker.k = @(d) var*(1 + sqrt(3) .* abs(d) ./ l) .* exp(-sqrt(3) .* abs(d) ./ l);
    ker.fam = 'matern32';
elseif nu == 2.5
    ker.k = @(d) var*(1 + sqrt(5).*abs(d)./l + 5.*abs(d).^2/(3*l^2)) .* exp(-sqrt(5) .* abs(d) / l);
    ker.fam = 'matern52';
end

% following taken from philip nufft_gps2d:  
scaling = (2*sqrt(pi))^dim * gamma(nu+dim/2) * (2*nu)^nu / (gamma(nu) * l^(2*nu));   % yuk; see https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
ker.khat = @(xid) var*scaling * (2*nu/l^2 + (4*pi^2) * xid.^2).^(-(nu + dim/2));
% xid means |xi|, ie the "distance" from origin in Fourier space.

% attributes to extract when using certain algorithms (e.g. ski)
ker.l = l;
ker.nu = nu;
ker.dim = dim;
