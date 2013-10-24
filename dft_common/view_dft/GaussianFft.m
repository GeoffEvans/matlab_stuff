% A unit gaussian aperture at the centre

FtOfGaussian = @(k)(gaussmf(k, [1, 0]) *  sqrt(2 * pi));
ViewFftWithFt( 256, 4, @(x) Gaussian(x), FtOfGaussian );
