% global matlab setup (will be automatically run if matlab started in this dir)

% point to library locations here, based on machine name
[status,host] = system('uname -n');
if strcmp(host(1:4),'ross')     % alex's setup
  addpath ~/numerics/finufft/matlab
  % for alex, MATLAB has to be started from conda env idp (py 3.7)
elseif strcmp(host(1:4),'C02G') % Manas' setup 
  addpath ~/git/finufft/matlab
else                       % default & all others (philip)
  addpath ~/nufft_gps/finufft/matlab;
end

% access all of local tree
h = fileparts(mfilename('fullpath'));
addpath(genpath(h))                        % gives access to all subdirs
rmpath(genpath(fullfile(h,'.git')))

% try prevent Py crashes at cost of performance...
% see: https://www.mathworks.com/matlabcentral/answers/486171-how-do-i-troubleshoot-a-matlab-crash-when-trying-to-use-the-python-interface?s_cid=pl_crsh_an
pyenv("ExecutionMode","OutOfProcess");

% add to python path for python calls
rel_path_to_ski = '/algs/SKI';
try
    path_to_ski = strcat(h, rel_path_to_ski);
    if count(py.sys.path, path_to_ski) == 0
        insert(py.sys.path, int32(0), path_to_ski);
    end
catch
    fprintf('Error in finding cython\n Skipping python path\n');
end

if(exist('algs/GPMLE/FLAM/startup.m'))
    run 'algs/GPMLE/FLAM/startup.m'
end

% terminal style
format long g
format compact

% restore pre-R2014b figure look:
set(groot,'DefaultFigureGraphicsSmoothing','off')   % kill slow antialiasing
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                                   .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))
set(groot,'defaultLegendAutoUpdate','off')
