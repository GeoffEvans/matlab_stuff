function [Y] = ViewFftWithFt( numberOfSamples, lengthOfAperture, ApertureFunction, FourierTransform )

intervalWidth = lengthOfAperture / numberOfSamples;
K = ((0:numberOfSamples-1) - numberOfSamples/2) * 2 * pi / lengthOfAperture;              % The wavevector samples
X = (0:numberOfSamples-1) * intervalWidth - lengthOfAperture/2;     % The position samples
Y = ApertureFunction(X);                                            % Image of the aperture
F = fft(SwapArrayHalves(Y));                                        % FFT of the aperture adjusted for phase
FN = real(F) * intervalWidth;                                       % Extract the amplitude and normalise 
FP = plot(K, SwapArrayHalves(FN));                             % Plot first the FFT
hold;
SP = plot(K, FourierTransform(K));                             % Compare to the actual FT
set(SP, 'color', 'red');
R = plot(K, (abs(K) * intervalWidth < 1) / 10 - 0.05);            % Approx is good for positive
set(R, 'color', 'green');
hold;
end

