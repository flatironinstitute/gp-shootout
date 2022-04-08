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

% add to python path for python calls
rel_path_to_ski = '/algs/SKI';
path_to_ski = strcat(h, rel_path_to_ski);
py_path = py.sys.path;
if count(py_path, path_to_ski) == 0
    insert(py_path, int32(0), path_to_ski);
end