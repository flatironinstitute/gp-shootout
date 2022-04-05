% global matlab setup

% terminal style
format long g
format compact

% point to library locations here, based on machine name
[status,host] = system('uname -n');
if strcmp(host(1:4),'ross')     % alex's setup
  addpath ~/numerics/finufft/matlab
else                       % default & all others (phillip)
  addpath ~/nufft_gps/finufft/matlab;
end

% access all of local tree
h = fileparts(mfilename('fullpath'));
addpath(genpath(h))                        % gives access to all subdirs
rmpath(genpath(fullfile(h,'.git')))
