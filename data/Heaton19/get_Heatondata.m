function [x, meas, dummy, xtrg, truetrg] = get_Heatondata(type)
% GET_HEATONDATA   return Heaton et al. '19 comparison data
%
% [x, meas, ~, xtrg, truetrg] = get_Heatondata(type) loads into RAM either
%  simulated (type='sim') or satellite measured (type='sat') 2D temperature
%  dataset. The sim data comes from a Matern-1/2 GP.
%  The coordinates are the raw data (LON,LAT), ie longitude in degrees and
%  latitude in degrees, without the aspect ratio correction that should be
%  applied to such a patch of the sphere; this is as in Heaton et al
%  publication. Data is mean surface (not air) temperature in degree Celcius.
%
%  See: data/Heaton19/README and make_Heaton_mat
%
% Outputs:
%  x    - 2*N array of coords of N data (training) points in R^2.
%  meas - length-N vector of temperature data (training).
%  [3rd output unused]
%  xtrg - 2*p array of coords of N target (test) points in R^2.
%  truetrg - length-p vector of "ground-truth" temperatures at test points.
%
% Without input or output arguments, does self-test.
if nargin==0 && nargout==0, test_get_Heatondata; return; end

if nargin==0, type='sim'; warning('using sim Heaton data by default'); end
h = fileparts(mfilename('fullpath'));      % dir this script is in
if strcmp(type, 'sim')
  load([h '/heaton19sim.mat'])
elseif strcmp(type, 'sat')
  load([h '/heaton19sat.mat'])
else
  error('unknown type!');
end
dummy = [];

%%%%%%%%
function test_get_Heatondata         % only visual for now (not unit test)
[x, meas, ~, xtrg, truetrg] = get_Heatondata('sat');
figure;
subplot(2,1,1); geoscatter(x(2,:), x(1,:), 1, meas); % geoscat expects LAT, LON
cb = colorbar; c = caxis; title('Heaton train');
subplot(2,1,2); geoscatter(xtrg(2,:), xtrg(1,:), 1, truetrg);    % "
caxis(c); cb = colorbar; title('Heaton test');
