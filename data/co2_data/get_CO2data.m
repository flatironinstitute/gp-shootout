function [x, meas] = get_CO2data(Nsub)
% GET_C02DATA   pull in co2 data set
%
% [x, meas] = get_CO2data() pulls in data from a C02 data set described
% here -- https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1419136
% (see Figure 4)
%
% Inputs:
%  Nsub - number of points to subsample from the ~1.5 million in data set
% Outputs:
%  x - 2 * N real array of latitude/longitude coordinates
%  meas - N * 1 real array of measured (observed) values
%
% Without arguments does self-test.
if nargin==0, test_get_CO2data; return; end
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


%%%%%%%%
function test_get_CO2data         % only visual for now (not unit test)
[x, meas] = get_CO2data(10000);
figure;
geoscatter(x(:,1), x(:,2), 1, meas);
cb = colorbar;

