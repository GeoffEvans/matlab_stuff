function [ k,ft ] = dft( func, xResolution, xStart, xNumSamples )

xLength = xResolution * xNumSamples;
kResolution = 2*pi/xNumSamples / xResolution;

x = xResolution * (0:xNumSamples-1);
xPeriodic = mod(x-xStart,xLength) + xStart;
ft = fft(func(xPeriodic));

kVals = kResolution * (0:xNumSamples-1);
kCutoff = kResolution * xNumSamples;
k = kVals - (kVals >= kCutoff/2) * kCutoff;

end

