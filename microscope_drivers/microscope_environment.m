wavelength = 800 * 1e-9;

optimalFrequency = 40 * 1e6;
acceptanceAngle = 4.35 * 1e-3; % rad

systemClockFreq = 240e6;
dataTimeInterval = 1 / 20e6;
zoomFactor = 1;

xyNumOfElems = 100;

dwellTimeMultiplier = 1;
xImageCentreNormalised = 0; 
yImageCentreNormalised = 0; 
zImageCentreNormalised = 0;

pairDeflectionRatio = 1;

[A, B, C, tRamp] = microscope_drivers_structural(...
    xyNumOfElems,...
    acceptanceAngle,... % mrad
    dwellTimeMultiplier,...
    wavelength,... % currently set in aod3d.m, can't imagine this ever changes
    zoomFactor,...
    optimalFrequency,...
    xImageCentreNormalised,...
    yImageCentreNormalised,...
    zImageCentreNormalised,...
    systemClockFreq,...
    pairDeflectionRatio);
    

