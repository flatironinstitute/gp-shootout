% print for tikzpicture
function print_tikz(xs, ys)
    n = numel(ys);
    for i=1:n
        fprintf('(%.3g, %.3g)\n', xs(i), ys(i));
    end
end