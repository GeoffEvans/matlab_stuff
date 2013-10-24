%Function length is 551cm with a sinc aperture of period 1/4pi cm focused at the centre

FtOfSinc = @(K) 2*pi*Boxcar(K);
ViewFftWithFt( 2048, 551, @(x) sinc(x/(2*pi)), FtOfSinc);
