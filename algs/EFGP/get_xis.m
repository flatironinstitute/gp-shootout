function [xis, h, mtot] = get_xis(ker, eps, L, opts)
% GET_XIS   Return 1D equispaced Fourier quadrature nodes for given tolerance
%
% xis = get_xis(ker, eps, L) returns an equispaced list of 1D frequency nodes
%  xis ("\xi's") adequate for Fourier quadrature in the EFGP method for the
%  kernel ker, to absolute tolerance eps, assuming that all coordinate
%  differences lie in [-L,L].
%
%  Inputs:
%   ker   - kernel struct with at least fields k, khat, and dim.
%   eps   - tolerance parameter, eg 1e-6
%   L     - max size of spatial domain in any coordinate, so that all
%           inter-point differences lie in [-L,L].
%  opts   - Optional arguments with default values in parenthesis
%              opts.use_integral (false), use the integral estimate
%               heuristic by philip
%              opts.l2scaled (false) whether to use l2 scaling of integral 
%               with new heuristics in determining h,m
%  Output:
%   xis   - row-vector of equispaced values to be used as quadrature nodes
%           (all weights in 1D are h=xis(2)-xis(1)).  Their number is odd.
%   h     - spacing between nodes
%   mtot  - total number of nodes in 1D (2*m+1 in paper)
%  Notes:
%  1) The spacing h is chosen using L, the real-space radial kernel function
%     ker.k, the spatial dimension ker.dim, and an aliasing error estimate.
%     Alex changed this to >Nyquist based on aliasing theory, 4/19/22.
%  2) The cutoff max(abs(xis)) is chosen via the relative decay of ker.khat.
%
% To test, see EFGP.
  
  use_integral = false;
  if(nargin == 3)
      opts = [];
  end
  if(isfield(opts,'use_integral'))
      use_integral = opts.use_integral;
  end
  dim = ker.dim;              % spatial dimension
  
  % spatial radial ker func
  k = ker.k;
  khat = ker.khat;
  eps_use = eps;
 
  
  if(use_integral)
  
      Ltime = getL(eps, k);       % find eps-support
      h = 1/(L+Ltime);            % xi node spacing so nearest aliased tail <= eps
      % (NB Nyquist was hnyq = 1/(2*Ltime), more pessimistic)

      % Fourier radial ker func

      %%%khat = @(r) ker.khat(r) / ker.khat(0);         % some old version w/o polar
      khat = @(r) abs(r^(dim-1)) * khat(r) / khat(0);   % polar factor & rel to 0
      Lfreq = getL(eps, khat);    % find eps-support

      hm = ceil(Lfreq/h);         % half number of nodes to cover [-Lfreq,Lfreq]
  else
      if(contains(ker.fam,'matern'))
          l = ker.l;
          nu = ker.nu;
          dim = ker.dim;
          eps_use = eps/ker.var;
          if(isfield(opts,'l2scaled')) 
             if(opts.l2scaled)
               % the following is the L2 norm of the kernel k (see p. 24 of the original arxiv paper) 
               rl2sq = (2*nu/pi/l^2)^(dim/2)*ker.khat(0)^2/2*gamma(dim/2+2*nu)/gamma(dim+2*nu)*2^(-dim/2);  % alex notes there's cancellation of 2^stuff here?
               eps_use = eps*sqrt(rl2sq);
             end
          end

          eps = eps_use;
          h = 1/(L+0.85*l/sqrt(ker.nu)*log(1/eps));   % heuristic \eqref{hheur}
          % note hm is "m" in the paper...
          hm = ceil(( pi^(nu+dim/2)*l^(2*nu) * eps/0.15 )^(-1/(2*nu+dim/2)) / h); % heuristic \eqref{mheur}
          
      elseif(strcmpi(ker.fam,'squared-exponential'))
          l = ker.l;
          dim = ker.dim;
          var = ker.var;
          eps_use = eps/var;
          if(isfield(opts,'l2scaled'))
             if(opts.l2scaled)
               rl2sq = ker.k(0)^2*(sqrt(pi)*l^2)^dim;
               eps_use = eps*sqrt(rl2sq);
             end
          end
          eps = eps_use;
          h = 1/(L+l*sqrt(2*log(4*dim*3^dim/eps)));
          hm = ceil(sqrt(log(dim*(4^(dim+1))/eps)/2)/pi/l/h); % again, "m"
      end
  end
  
  
  xis = (-hm:hm)*h;           % use exactly h, so can get bit of spillover
  mtot = numel(xis);          % 2m+1
end
