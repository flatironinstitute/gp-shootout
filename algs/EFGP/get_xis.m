function [xis, h, m] = get_xis(ker, eps, L)
% GET_XIS   Return 1D equispaced Fourier quadrature nodes for given tolerance
%
% xis = get_xis(ker, eps, L) returns an equispaced list of 1D frequency nodes
%  xis ("\xi's") adequate for Fourier quadrature in the EFGP method for the
%  kernel ker, to absolute tolerance eps, assuming that all coordinate
%  differences lie in [-L,L].
%
% [xis h m] = get_xis(ker, eps, L) also returns spacing h and number of nodes m.
%
%  Inputs:
%   ker   - kernel struct with at least fields k, khat, and dim.
%   eps   - tolerance parameter, eg 1e-6
%   L     - max size of spatial domain in any coordinate, so that all
%           inter-point differences lie in [-L,L].
%  Output:
%   xis   - row-vector of equispaced values to be used as quadrature nodes
%           (all weights in 1D are h=xis(2)-xis(1)).  Their number is odd.
%   h     - spacing
%   m     - number of nodes
%  
%  Notes:
%  1) The spacing h is chosen using L, the real-space radial kernel function
%     ker.k, the spatial dimension ker.dim, and an aliasing error estimate.
%     Alex changed this to >Nyquist based on aliasing theory, 4/19/22.
%  2) The cutoff max(abs(xis)) is chosen via the relative decay of ker.khat.
%
% Testing: for now, run EFGP.
  
  dim = ker.dim;              % spatial dimension
  
  % spatial radial ker func
  k = ker.k;                  % care about absolute vals of kernel
  Ltime = getL(eps, k);       % find eps-support
  h = 1/(L+Ltime);            % xi node spacing so nearest aliased tail <= eps
  % (NB Nyquist was hnyq = 1/(2*max(L,Ltime)), more pessimistic)
  
  % Fourier radial ker func
  khat = ker.khat;
  %%%khat = @(r) ker.khat(r) / ker.khat(0);         % some old version w/o polar
  khat = @(r) abs(r^(dim-1)) * khat(r) / khat(0);   % polar factor & rel to 0
  Lfreq = getL(eps, khat);    % find eps-support
  
  hm = ceil(Lfreq/h);         % half number of nodes to cover [-Lfreq,Lfreq]
  xis = (-hm:hm)*h;           % use exactly h, so may have bit of spillover
  m = numel(xis);
end
