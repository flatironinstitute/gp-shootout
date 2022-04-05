function r = rms(x)
% RMS   root mean square of all elements of an array.
r = norm(x(:))/sqrt(numel(x));
% to test: rms(rand(1e4,1)) should be 1 +- 0.01-ish
