%Function length is 15cm with a 1cm boxcar aperture at the centre

FtOfBox = @(k)(sinc(k/(2*pi)));
ViewFftWithFt( 4096, 15, @(x) Boxcar(x), FtOfBox);
