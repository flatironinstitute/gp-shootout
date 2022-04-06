function L = getL(eps, kern)
%
% find the value of kern that's eps relative to its maximum. we assume that
% kern, has a maximum at 0 and decays monotonically. this is done using
% bisection. 
%
% Inputs:
% eps - value at which to truncate
% kern - kernel 
%
% Outputs:
% L - value at which kern(L) = eps
%

    % make sure starting upper bound is large enough
    a = 0;
    b = 1000;
    nmax = 10;
    for i=1:nmax
        if kern(b) > eps
            b = 2*b;
        else
            break
        end
    end

    % start bisection
    nmax = 200;
    for i=1:nmax
        mid = (a + b)/2;
        fmid = kern(mid);
        if fmid > eps
            a = mid;
        else
            b = mid;
        end
    end
    L = mid;
end