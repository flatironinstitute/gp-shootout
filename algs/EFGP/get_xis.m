function [xis] = get_xis(ker, eps, tmax)
% for a given kernel, find the equispaced fourier discretization of that
% kernel. do this by truncating the kernel in time domain, then using the
% nyquist spacing in frequency domain, truncating the expansion when the
% magnitude of the fourier transform is eps. 

    % find support of function in time domain
    Ltime = getL(eps, ker.k);
    Ltime = max(tmax, Ltime);
    
    % nyquist
    hnyq = 1/(2*Ltime);
    
    % find numerical bandlimit via bisection
    Lfreq = getL(eps, ker.khat);
    
    % number of nodes for nyquist, should be odd
    m = 2*Lfreq/hnyq + 1;
    m = 2*ceil(m/2) + 1;
    xis = linspace(-Lfreq, Lfreq, m);
end