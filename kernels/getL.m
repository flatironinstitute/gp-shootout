function L = getL(eps, kern)
% GETL   finds where monotonically descreasing function reaches given value
%
% L = getL(eps, kern) finds L where kern(L)=eps, assuming that kern is a handle
% to a monotonically decreasing scalar function.
% The algorithm is geometric expansion followed by bisection search.
%
% Inputs:
% eps - absolute value at which to truncate
% kern - function handle (eg, kernel or its Fourier transform)
%
% Outputs:
% L - approximate value at which kern(L) = eps


%  eps = eps*kern(0);   % absolute value to match - AHB tried
  
  % make sure starting upper bound is large enough
  a = 0;
  b = 1000;
  nmax = 10;
  for i=1:nmax
    if kern(b) > eps
      b = 2*b;
    else
      break
    end
  end
  
  % start bisection
  nmax = 200;
  for i=1:nmax
    mid = (a + b)/2;
    fmid = kern(mid);
    if fmid > eps
      a = mid;
    else
      b = mid;
    end
  end
  L = mid;
end
