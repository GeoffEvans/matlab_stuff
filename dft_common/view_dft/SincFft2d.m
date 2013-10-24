% A sinc aperture of period 1/4pi cm along x and y directions both focused
% at the origin

sinc2d = @(x, y) sinc(x/(2*pi)) .* sinc(y/(2*pi));
ViewFft2d( 1024, 851, sinc2d);
