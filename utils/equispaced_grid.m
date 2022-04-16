function xs = equispaced_grid(dim, n)
% construct dim x n grid of points on [0, 1]^n
    if dim == 1
        xs = linspace(0, 1, n);
    end

    if dim == 2
        x = linspace(0, 1, n);
        [xx, yy] = meshgrid(x);
        xs = [reshape(xx, 1, []), reshape(yy, 1, [])];
    end

    if dim == 3
        x = linspace(0, 1, n);
        [xx, yy, zz] = meshgrid(x);
        xs = [reshape(xx, 1, []), reshape(yy, 1, []), reshape(zz, 1, [])];
    end

end

