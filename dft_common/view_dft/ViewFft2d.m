function [Y] = ViewFft2d( numberOfSamples, lengthOfAperture, ApertureFunction2d )

intervalWidth = lengthOfAperture / numberOfSamples;
K = ((0:numberOfSamples-1) - numberOfSamples/2) * 2 * pi / lengthOfAperture;	% The wavevector sample basis
[Kx, Ky] = meshgrid(K);                                                         
R = (0:numberOfSamples-1) * intervalWidth - lengthOfAperture/2;     % The position sample basis
[X, Y] = meshgrid(R);
A = ApertureFunction2d(X, Y);                                       % 2D Image of the aperture
F = fft2(SwapMatrixHalves(A));                                      % 2D FFT of the aperture adjusted for phase
FN = SwapMatrixHalves(real(F));                                      % Extract the amplitude and normalise 
goodApprox = (abs(Kx) * intervalWidth < 1) & (abs(Ky) * intervalWidth < 1);   % Approx is good for zero
GA = surf(Kx, Ky, goodApprox * F(1,1));                                         % Plot region of good approx
alpha(GA, 0.2)
set(GA,'LineStyle','none');
hold on;
FP = surf(Kx, Ky, FN);                                              % Plot FFT
alpha(FP, 0.7)
set(FP,'LineStyle','none');
hold off;
end

