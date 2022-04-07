% matlab by default uses version 2.7, but the Python code called was
% written in version 3.7. users may need to replace the path with
% the location of python3 on their machinies

% if we haven't successfully set the version to 3.7, do that now
pe = pyenv;
if ~strcmp(pe.Version, '3.7')
    path_to_py3 = '/usr/local/bin/python3';
    pyenv('Version',path_to_py3);
end

% reload python module in case changes were made to python file. this is
% useful when making changes to the python file so that you don't have to
% restart matlab each time
%clear classes
%m = py.importlib.import_module('ski1d');
%py.importlib.reload(m);


% data
dim = 1;
N = 10;
x = linspace(0, 1, N)';
x = rand(N, dim);
x(1) = 0;
x(N) = 1;
x = sort(x);
y = cos(x) + rand(N, 1);

% test points
ntest = 100;
%testx = linspace(0, 1, ntest)';
testx = sort(rand(ntest, dim));
testx(1) = x(1);
testx(ntest) = x(N);
grid_size = 100;
double_prec = 1;

l = 0.1;
sigma2 = 1.0;
% kern_family must be one of 'matern12', 'matern32', 'squared-exponential'
kern_family = 'matern12';
kern_family = 'squared-exponential';
xpy = py.numpy.array(x);
ypy = py.numpy.array(y);
testxpy = py.numpy.array(testx);
out1 = py.ski.ski_mat(xpy, ypy, testxpy, grid_size, sigma2, kern_family, l);

% unpack ski output
yhat = double(out1{1})';
time_info = out1{2};

% compare to naive approach
ker = SE_ker(1, l);
[yhat2, ytrg, info] = naive_gp(x', y, sigma2, ker, testx', []);

disp(yhat - ytrg.mean);

% scatter(testx, yhat);
% hold on
% scatter(testx, ytrg.mean);
% scatter(x, y);
% hold off




%%% 2d example
% data
dim = 2;
N = 10;
x = rand(N, dim);
x(1, :) = [0, 0];
x(N, :) = [1, 1];
y = cos(x(:, 1) + x(:, 2)) + rand(N, 1);

% test points
ntest = 10;
testx = rand(ntest, dim);
%testx(1) = x(1);
%testx(ntest) = x(N);
grid_size = 100;
double_prec = 1;

l = 0.1;
sigma2 = 1.0;
% kern_family must be one of 'matern12', 'matern32', 'squared-exponential'
kern_family = 'matern12';
kern_family = 'squared-exponential';
xpy = py.numpy.array(x);
ypy = py.numpy.array(y);
testxpy = py.numpy.array(testx);
out1 = py.ski.ski_mat(xpy, ypy, testxpy, grid_size, sigma2, kern_family, l);

% unpack ski output
yhat = double(out1{1})';
time_info = out1{2};

% compare to naive approach
ker = SE_ker(1, l);
[yhat2, ytrg, info] = naive_gp(x', y, sigma2, ker, testx', []);

disp(yhat - ytrg.mean);
