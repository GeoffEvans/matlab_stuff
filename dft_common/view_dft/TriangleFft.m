%Function length is 15cm with a 1cm triangle aperture at the centre

FtOfTriangle = @(k)(sinc(k/(4*pi)).^2 ./ 2);
ViewFftWithFt( 4096, 15, @(x) Triangle(x), FtOfTriangle);
