function [ k,ft ] = dft( func, xResolution, xStart, xEnd )

numSamples = 2^ceil(log2((xEnd - xStart)/xResolution)); % round to even num
xLength = xResolution * numSamples;
xEnd = xStart + xLength;
kLength = 2 * pi * numSamples / xLength;
kResolution = 2*pi/numSamples / xResolution;

x = xResolution * (0:numSamples-1);
xPeriodic = mod(x-xStart,xLength) + xStart;
ft = fft(func(xPeriodic));

kVals = kResolution * (0:numSamples-1);
kCutoff = kResolution * numSamples;
k = kVals - (kVals >= kCutoff/2) * kCutoff;

end

