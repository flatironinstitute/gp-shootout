% run all tests of pure-matlab things.
% See README.md
% (Also see: test_all_nonmatlab for trickier cross-language interfaces).

test_kers
densekermat
get_randdata
naive_gp
try
    EFGP
catch
    fprintf('EFGP not installed successfully\n');
end

try
    FLAMGP
catch
    fprintf('FLAGMP not installed successfully\n');
end
