function s=scifmt(x,d)
% SCIFMT  Convert real number x into exponent format string with d sig figs

% Barnett 10/13/22
  if nargin==0, test_scifmt; return; end
  if nargin<2, d=2;  end   % default
  
  exp = log10(x);      % hack the eps latex string  ... same for all N
  iexp = floor(exp);
  if x==10^iexp         % only if exact to machine prec --- watch out
    s = sprintf('10^{%d}',iexp);
  else
    fmt = sprintf('%%.%df',d-1);       % meta format
    fmt = [fmt '\\times 10^{%d}'];     % escaped \
    s = sprintf(fmt, x/10^iexp, iexp);
  end

  %%%%%
  function test_scifmt
    scifmt(0.01)
    scifmt(0.02)
    scifmt(pi*1e-6)
    scifmt(pi*1e-6,3)

    