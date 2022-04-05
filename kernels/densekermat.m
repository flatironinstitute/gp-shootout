function K = densekermat(k,x,y)
% DENSEKERMAT   Fill dense matrix using a isotropic kernel function
%
% densekermat(k,x) where x is d*m matrix of d-dimensional coords of M points,
%  returns M*M symmetric dense matrix with ij'th element k(||x_i - x_j||).
%
% densekermat(k,x,y) where additionally y is a d*N matrix of source points
%  returns the M*N dense matrix with ij'th element k(||x_i - y_j||). Here x
%  is interpreted as target points, so that the matrix maps R^N to R^M.
%
% Self-test done if called without args
if nargin==0, test_densekermat; return; end

if nargin<3, y=x; end
dim = size(x,1);
if size(y,1)~=dim, error('y must have same # dims as x!'); end
dd = abs(y(1,:)-x(1,:)');     % matrix of abs diffs, note broadcast op
if dim==1
  K = k(dd);
else
  dd = dd.^2;                 % now dd = distance^2 matrix; add dim-by-dim
  for d=2:dim
    dd = dd + (y(d,:)-x(d,:)').^2;
  end
  K = k(sqrt(dd));
end


%%%%%%%
function test_densekermat     % size test only, not vals!
N = 1e2; % src
M = 30;  % trg
f = @(d) sin(abs(d));
for dim=1:2
  x = rand(dim,N);
  A = densekermat(f,x);
  B = densekermat(f,rand(dim,M),x);
  assert(size(A,1)==N && size(A,2)==N && norm(A-A')==0.0)
  assert(size(B,1)==M && size(B,2)==N)
end
disp('densekermat test passed')
