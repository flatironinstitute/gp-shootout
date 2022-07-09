% test all things which require compilation or inter-language pain
try
    SKI
catch
    fprintf('SKI not installed successfully\n');
end

try
    RLCM
catch
    fprintf('RLCM not installed successfully\n');
end
