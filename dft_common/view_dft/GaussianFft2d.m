% A 2d unit gaussian aperture at the centre

gauss2d = @(x, y) gaussmf(x, [1 0]) .* gaussmf(y, [1 0]);
ViewFft2d( 1024, 351, gauss2d);
