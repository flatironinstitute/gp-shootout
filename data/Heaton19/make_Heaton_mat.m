% script to examine and select from CSV files from Heaton et al. 2019
% comparison 2D satellite data, and write out to our MAT files.
% Optionally rescales to [0,1]^2 and fixes aspect ratio, extracts test
% targets and true values therein.
% See README for how we extracted .csv files from .RData (R format) files.
% Barnett 4/11/22. rescaling switched off 4/19/22.
%
% Notes:
% 1) it's clear that Heaton et al. used raw (LON,LAT) coords to construct
%    simulated data using an isotropic Matern kernel:
%    https://github.com/finnlindgren/heatoncomparison/blob/master/Code/FormatData/FormatData.R
%    However this is *not* locally isotropic on the Earth surface!
%    They omitted the sin(theta) factor, where theta = (90-LAT)*pi/180 is
%    radians from the north pole.
%    How could this simple thing be missed? It is intentional?
% 2) LST (land surface temperatures) are in deg Celcius, see:
%   https://ladsweb.modaps.eosdis.nasa.gov/filespec/MODIS/6/MOD11_L2
%   showing scale_factor = 0.02
%   used in line 25 of
%   https://github.com/finnlindgren/heatoncomparison/blob/master/Code/FormatData/FormatData.R
% 3) LST measured as in
%   https://lpdaac.usgs.gov/documents/119/MOD11_ATBD.pdf
%   and claims 0.5 to 1.0 K accuracy (hence also Celcius).

clear; verb = 1;      % verbosity (0=no figs, 1=figs)

% location of your clone of heatoncomparison repo...
p = '/home/alex/numerics/nufft/GP/heatoncomparison/Data/';
h = fileparts(mfilename('fullpath'));      % dir this script is in

ar=readmatrix([p 'AllSatelliteTemps.csv'],'numheaderlines',1);   % "r" = real
% cols: index longitude(deg) latitude(deg) alltemps traintemps.
% alltemps are Nan where no meas. traintemps are Nan where instead a test pt
% scalar temperature data appears to be in degrees F.
r=readmatrix([p 'SatelliteTemps.csv'],'numheaderlines',1);
as=readmatrix([p 'AllSimulatedTemps.csv'],'numheaderlines',1);   % "s" = sim
s=readmatrix([p 'SimulatedTemps.csv'],'numheaderlines',1);
% note: didn't need to load r and s files, since subset of ar and as :)

if verb>1, figure;  % show raw data in degree (long,lat) units...
  subplot(2,2,1); scatter(ar(:,2),ar(:,3),1,ar(:,5)); colorbar
  title(sprintf('ar (LON,LAT): all sat, N=%d',size(ar,1)));
  subplot(2,2,2); scatter(r(:,2),r(:,3),1,r(:,4)); colorbar
  title(sprintf('r (LON,LAT): sat, N=%d',size(r,1)));
  subplot(2,2,3); scatter(as(:,2),as(:,3),1,as(:,5)); colorbar
  title(sprintf('as (LON,LAT): all sim, N=%d',size(as,1)));
  subplot(2,2,4); scatter(s(:,2),s(:,3),1,s(:,4)); colorbar
  title(sprintf('s (LON,LAT): sim, N=%d',size(s,1)));
end

rescale = false;               % false: leave original x coords & aspect ratio.
correctaspect = false;         % Heaton uses false; we believe true is better

lon0 = 0.0; lat0 = 0.0; lonsc = 1.0; latsc = 1.0; sc = 1.0; sth = 1.0; % default
if rescale
  lon0 = min(ar(:,2));             % get bounding box around the data points
  lonsc = max(ar(:,2))-min(ar(:,2));
  lat0 = min(ar(:,3));
  latsc = max(ar(:,3))-min(ar(:,3));
  % project roughly to flat projection w/ correct aspect ratio...
  if correctaspect
    latcen = lat0+latsc/2;           % latitude of center of box
    sth = cos(2*pi*latcen/360);      % "sin(theta)" metric factor (lon vs lat)
  end
  sc = max(lonsc*sth,latsc);    % conversion from lat degrees so x in [0,1]^2
end
lonfac = sth/sc;              % factor we multiplied LON by to get x(1,:)
latfac = 1/sc;                % factor we multiplied LAT by to get x(2,:)
xa(1,:) = (ar(:,2)'-lon0) * lonfac;    % all x-coords, possibly rescaled
xa(2,:) = (ar(:,3)'-lat0) * latfac;    %   " y
REkm = 6371;                     % typ radius of Earth, kilometer (1e-3 relerr)
kmdeg = REkm*pi/180;             % km per latitude degree
fprintf('LON*LAT box size %.3g*%.3g deg (%.3g*%.3g km)\n',lonsc,latsc,kmdeg*lonsc*sth,kmdeg*latsc)
fprintf('pixel sizes in km: %.4g*%.4g (2nd only valid on actual Earth sphere if correctaspect=true)\n',diff(xa(1,1:2))*sc*kmdeg, diff(xa(2,[501 1]))*sc*kmdeg)
if verb>1, figure; plot(xa(1,:),xa(2,:),'.'); title('all nodes xa'); end

% select training and test targets in more sensible format... real satellite
jtrain = find(~isnan(ar(:,4)));                    % "training" indices
jtest = find(isnan(ar(:,4)) & ~isnan(ar(:,5)));   % test targets indices
xtrg = xa(:,jtest);
x = xa(:,jtrain);         % "x" our name for training pts
meas = ar(jtrain,5);      % satellite temperature data (Celcius)
truetrg = ar(jtest,5);    % "
save([h '/heaton19sat.mat'],'x','meas','xtrg','truetrg','lonfac','latfac');
% simulated (but note simulation used kernel incorrectly isotropic in LON,LAT)
jtrain = find(~isnan(as(:,4)));                    % "training" indices
jtest = find(isnan(as(:,4)) & ~isnan(as(:,5)));   % test targets indices
xtrg = xa(:,jtest);
x = xa(:,jtrain);         % "x" our name for training pts
meas = as(jtrain,5);      % simulated temperature data (Celcius). Same x, xtrg
truetrg = as(jtest,5);    % "
save([h '/heaton19sim.mat'],'x','meas','xtrg','truetrg','lonfac','latfac');

if verb, figure;  % show processed data
  load([h '/heaton19sat.mat'])
  subplot(2,2,1); scatter(x(1,:),x(2,:),1,meas); colorbar; axis equal tight;
  v = axis; c = caxis;
  title(sprintf('make\\_Heaton\\_mat: sat train N=%d',numel(meas)));
  subplot(2,2,2); scatter(xtrg(1,:),xtrg(2,:),1,truetrg); axis(v); caxis(c);
  colorbar; axis equal tight;
  title(sprintf('make\\_Heaton\\_mat: sat test Ntrg=%d',numel(truetrg)));
  load([h '/heaton19sim.mat'])
  subplot(2,2,3); scatter(x(1,:),x(2,:),1,meas); colorbar; axis equal tight;
  v = axis; c = caxis;
  title(sprintf('make\\_Heaton\\_mat: sim train N=%d',numel(meas)));
  subplot(2,2,4); scatter(xtrg(1,:),xtrg(2,:),1,truetrg); axis(v); caxis(c);
  colorbar; axis equal tight;
  title(sprintf('make\\_Heaton\\_mat: sim test Ntrg=%d',numel(truetrg)));
end
