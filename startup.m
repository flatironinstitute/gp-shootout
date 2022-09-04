% global matlab setup (will be automatically run if matlab started in this dir)

if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

iefgp = 0;
iflam = 0;
irlcm = 0;
igpytorch = 0;

str = which('finufft');   % checks the MEX file, not the matlab wrappers
if(~isempty(str))
   iefgp = 1;
end

% % point to library locations here, based on machine name
% [status,host] = system('uname -n');
% if strcmp(host(1:4),'ross')     % alex's setup
%   addpath ~/numerics/finufft/matlab
%   % for alex, MATLAB has to be started from conda env idp (py 3.7)
% elseif strcmp(host(1:4),'C02G') % Manas' setup 
%   addpath ~/git/finufft/matlab
% else                       % default & all others (philip)
%   addpath ~/nufft_gps/finufft/matlab;
% end

% access all of local tree
% try prevent Py crashes at cost of performance...
% see: https://www.mathworks.com/matlabcentral/answers/486171-how-do-i-troubleshoot-a-matlab-crash-when-trying-to-use-the-python-interface?s_cid=pl_crsh_an


% add to python path for python calls
rel_path_to_ski = '/algs/SKI';
try
    pyenv("ExecutionMode","OutOfProcess");
    path_to_ski = strcat(h, rel_path_to_ski);
    if count(py.sys.path, path_to_ski) == 0
        insert(py.sys.path, int32(0), path_to_ski);
    end
catch
    fprintf('Error in finding cython\n Skipping python path\n');
end



if(exist('algs/GPMLE/FLAM/startup.m','file'))
    run 'algs/GPMLE/FLAM/startup.m'
    iflam = 1;
else
    str = which('rskelf');
    if(~isempty(str))
        iflam = 1;
    end
end

igpp = 0;
if(exist('algs/RLCM','dir'))
    
    
    for j=12:-1:6
        struse = ['which g++-' num2str(j)];
        [status,cmdout] = system(struse);
        if(~status)
            fprintf('Found gnu g++ compiler version = %d\n',j);
            igpp = 1;
            break
        end
        irlcm = 1;
    end
    [status1,cmdout1] = system('which make');
    if(status1)
        igpp = 0;
    end
end

gpp = ['g++-' num2str(j)];
if(igpp && (ismac || isunix))
    path1 = getenv('PATH');
    cmdout2 = extractBefore(cmdout,gpp);
    path1 = [path1 cmdout2];
    setenv('PATH',path1);
    if(ismac)
        cd 'algs/RLCM';
        system(['./buildit_mac.sh ' gpp]);
        cd ../../;
    else
        cd 'algs/RLCM';
        system(['./buildit_linux.sh ' gpp]);
        cd ../../;
    end
end
    

h = fileparts(mfilename('fullpath'));
addpath(genpath(h))                        % gives access to all subdirs
rmpath(genpath(fullfile(h,'.git')))

% terminal style
format long g
format compact

% restore pre-R2014b figure look:
set(groot,'DefaultFigureGraphicsSmoothing','off')   % kill slow antialiasing
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                                   .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))
set(groot,'defaultLegendAutoUpdate','off')

fprintf('\n\n\n\n');
fprintf('---------------------------------------------\n');
if(iefgp)
    fprintf('EFGP: Installation successful\n');
else
    fprintf('EFGP: installation failure\n Finufft not found in path\n'); 
end
fprintf('---------------------------------------------\n');
if(iflam)
    fprintf('FLAM: Installation successful\n');
else
    fprintf('FLAM: installation failure\n Submodule not included and FLAM not found in path\n');
end
fprintf('---------------------------------------------\n');

if(irlcm)
    fprintf('RLCM: Installation successful\n');
else
    fprintf('RLCM: installation failure\n Submodule not included or g++ compiler or make not found\n');
end
fprintf('---------------------------------------------\n');

