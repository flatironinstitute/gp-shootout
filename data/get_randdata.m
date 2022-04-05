function [x, meas, truemeas] = get_randdata(dim, N, f, sigma)
% GET_RANDDATA   simple iid noisy random data in [0,1]^d from a function
%
% [x, meas, truemeas] = get_randdata(dim, N, f, sigma) generates a simulated
%  simple iid random dataset to test GP regression on. If dim=1 then the
%  data is sorted in x for easy plotting (we assume the algorithm is permutation
%  invariant).
%
% Inputs:
%   dim - dimension
%   N - number of points
%   f - function handle (vectorized) mapping all points to noise-free values
%   sigma - additive (measurement) noise level
% Outputs:
%  x - dim*N real array of point coordinates in [0,1]^dim
%  meas - N*1 real array of measured (observed) values
%  truemeas - N*1 real array of underlying values before noise added
%
% Without arguments does self-test.
if nargin==0, test_get_randdata; return; end
x = rand(dim,N);     % in [0,1]^dim
if dim==1, x = sort(x); end     % just to help plotting
truemeas = f(x);
meas = truemeas(:) + sigma*randn(N,1);   % noisy data, col vec


%%%%%%%%
function test_get_randdata          % only visual for now (not unit test)
N=1e4;
dim=2;
f = @(x) cos(6*(x(1,:).^2 + x(2,:).^2));      % squared-dist func in 2D
sigma = 0.5;
[x, meas, truemeas] = get_randdata(dim, N, f, sigma);
figure;
subplot(1,2,1); scatter(x(1,:),x(2,:),[],truemeas,'filled');
caxis([-1 1]); axis equal tight
subplot(1,2,2); scatter(x(1,:),x(2,:),[],meas,'filled');
caxis([-1 1]); axis equal tight
title('get\_randdata test');
