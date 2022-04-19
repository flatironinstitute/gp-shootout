function [x, meas] = get_CO2data(Nsub)
% GET_C02DATA   pull in CO_2 2D data set.
%
% [x, meas] = get_CO2data pulls in data from a C02 data set described
% here -- https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1419136
% (see Figure 4).
%
% [x, meas] = get_CO2data(Nsub) restricts to a random subset of size Nsub.
%
% Inputs:
%  Nsub - [optional] number of points to subsample from the ~1.5 million in
%         original data set
% Outputs:
%  x    - 2*N real array of latitude/longitude coordinates
%  meas - N*1 real array of measured (observed) values, CO_2 conc in ppm.
%
% Without input or output arguments, does self-test.
if nargin==0 && nargout==0, test_get_CO2data; return; end

fileID = fopen('co2_meas.bin');
meas = fread(fileID, 'double');
fclose(fileID);
N = numel(meas);

fileID = fopen('co2_xs.bin');
x = fread(fileID, [N 2], 'double');
fclose(fileID);

% uniformly subsample data 
if Nsub < N
    inds = randperm(N, Nsub);
    x = x(inds, :);
    meas = meas(inds);
end

x = x';           % to match doc
x = x([2 1],:);   % swap rows to (LON,LAT) format as with Heaton


%%%%%%%%
function test_get_CO2data         % only visual for now (not unit test)
[x, meas] = get_CO2data(100000);
figure;
geoscatter(x(2,:), x(1,:), 1, meas);      % geoscat expects LAT, LON
cb = colorbar; title('CO2 data subset (ppm)');
