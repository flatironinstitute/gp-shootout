function ker = SE_ker(dim,l)
% SE_KER   Squared-exponential kernel func and Fourier transform, general dim
%
% ker = SE_ker(dim,l) returns ker.k a function handle taking an array of
%  distances to kernel values, and ker.khat, the function handle to its
%  Fourier transform (mapping wavevector magnitudes to values).
%  Latter uses Phillip's convention. For squared-exponential kernel (Gaussian),
%  l is the distance scale.
%
% To test see: TEST_KERS

if dim<1, error('dim must be at least 1!'); end
if l<=0, error('l lengthscale must be positive!'); end
ker.k = @(d) exp( (-0.5/l^2) * d.^2);                                % d = |x|
ker.khat = @(xid) (2*pi*l^2)^(dim/2) * exp( (-2*pi^2*l^2) * xid.^2); % |xi|
