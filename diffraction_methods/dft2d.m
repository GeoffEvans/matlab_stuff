function [ kx,ky,ft ] = dft2d( func, xResolution, xStart, xNumSamples, yResolution, yStart, yNumSamples )

xLength = xResolution * xNumSamples;
kxResolution = 2*pi/xNumSamples / xResolution;

yLength = yResolution * yNumSamples;
kyResolution = 2*pi/yNumSamples / yResolution;

x = xResolution * (0:xNumSamples-1);
xPeriodic = mod(x-xStart,xLength) + xStart;
y = yResolution * (0:yNumSamples-1);
yPeriodic = mod(y-yStart,yLength) + yStart;
[xGrid, yGrid] = meshgrid(xPeriodic,yPeriodic);
ft = fft2(func(xGrid,yGrid));

kxVals = kxResolution * (0:xNumSamples-1);
kxCutoff = kxResolution * xNumSamples;
kx = kxVals - (kxVals >= kxCutoff/2) * kxCutoff;

kyVals = kyResolution * (0:yNumSamples-1);
kyCutoff = kyResolution * yNumSamples;
ky = kyVals - (kyVals >= kyCutoff/2) * kyCutoff;

[kx,ky] = meshgrid(kx,ky);

end

