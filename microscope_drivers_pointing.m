function aod_dev_script()

wavelength = 800 * 1e-9;
baseFrequency = 40 * 1e6;
aodMode = -1;
acceptanceAngle = 4.35 * 1e-3; % rad
systemClockFreq = 240e6;
dataTimeInterval = 1 / 20e6;
dwellTimeMultiplier = 1;
dwellTime = dwellTimeMultiplier * dataTimeInterval;

zoomInteger = 1;
xyElements = 100;
zElements = 1;

subdividedVox = xyElements * zoomInteger;
imageCentres = [0 0 0];

[opticalDist1,opticalDist2,opticalDist3] = SpecifyOpticalDistances();

V = 612.9;
aodAperture = 15e-3;
aodFillTime = dataTimeInterval * ceil( aodAperture / V / dataTimeInterval ); % rounded up to a multiple of dataTimeInterval
halfFreqShiftMax = acceptanceAngle * V / wavelength;
chirpMax = 2 * halfFreqShiftMax / aodFillTime;
opticalDistMinimum4 = V^2 / (2 * wavelength * chirpMax);

numOfImageCentresPerCycle = size(imageCentres,1);
numberOfCycles=round(1);

PMCtime=dataTimeInterval*ceil(1e-3/dataTimeInterval) % Pointing mode cycle time (= integer*dataclock s )
Scantime=numberOfCycles*PMCtime

[xGrid,yGrid,zGrid] = meshgrid(1:xyElements,1:xyElements,1:zElements);

for u=1:numOfImageCentresPerCycle;
    Mxgp(u) = floor( Xp(u) * subdividedVox / 2) + ceil(subdividedVox / 2); % zoomed voxel identifier
    Mygp(u) = floor( Yp(u) * subdividedVox / 2) + ceil(subdividedVox / 2);
    Mzgp(u) = floor( Zp(u) * subdividedVox / 2) + ceil(subdividedVox / 2);
    
    Mxp{u} = Mxgp(u) + xGrid - ceil(xyElements/2);
    Myp{u} = Mygp(u) + yGrid - ceil(xyElements/2);
    Mzp{u} = Mzgp(u) + zGrid - ceil(zElements/2);
    
    Xnp{u} = (2*Mxp{u} - subdividedVox + 1)/(subdividedVox - 1);
    Ynp{u} = (2*Myp{u} - subdividedVox + 1)/(subdividedVox - 1);
    Znp{u} = (2*Mzp{u} - subdividedVox + 1)/(subdividedVox - 1);
    
    imageCentres(u,3) = imageCentres(u,3) + (imageCentres(u,3) == 0) * 1e-6;   % to prevent divide by zero
    
    opticalDist4 = opticalDistMinimum4 / imageCentres(u,3); % ???
    
    chirp(1) = aodMode * (V^2/wavelength) / ( 2 * (opticalDist1+opticalDist2+opticalDist3+opticalDist4) );
    chirp(2) = aodMode * (V^2/wavelength) / ( 2 * (opticalDist2+opticalDist3+opticalDist4) );
    chirp(3) = aodMode * (V^2/wavelength) / ( 2 * (opticalDist3+opticalDist4) );
    chirp(4) = aodMode * (V^2/wavelength) / ( 2 * opticalDist4 );
    
    timeZeroFreqOffset(1)=-aodMode*2*(opticalDist3+opticalDist4)*opticalDist4/((opticalDist4+opticalDist3)*(opticalDist1+opticalDist2+2*(opticalDist3+opticalDist4)))*(Xnp(u)+skewZX.*Znp(u))*halfFreqShiftMax;%minus
    timeZeroFreqOffset(2)=aodMode*2*opticalDist4/(2*opticalDist4+opticalDist3+opticalDist2)*(Ynp(u)+skewZY.*Znp(u))*halfFreqShiftMax;% AOD Y1 correct I think 150408
    timeZeroFreqOffset(3)=aodMode* opticalDist4/(opticalDist4+opticalDist3)*(Xnp(u)+skewZX.*Znp(u))*halfFreqShiftMax;% AOD X2 correct I think 150408%minus
    timeZeroFreqOffset(4)=-aodMode* (Ynp(u)+skewZY.*Znp(u))*halfFreqShiftMax;% AOD Y2 correct I think 150408
    
    centreFreqOffset = timeZeroFreqOffset - chirp * aodFillTime/2;
    
    A(u,:) = ScaleAndRound(baseFrequency + centreFreqOffset(1:4), 2^32 / systemClockFreq);
    B(u,:) = ScaleAndRound(chirp, 2^32 / systemClockFreq / 8); %CHANGED to scale B as for scanning mode
    C(u,:) = zero(numOfImageCentresPerCycle,4);
    
    scanTime = aodFillTime + dwellTime; % the scan time for each point
    ticksPerRamp(u,:)=  repmat( ScaleAndRound(scanTime, systemClockFreq), 1, 4);
    freqStart(u,:) = baseFrequency + timeZeroFreqOffset - chirp*scanTime/2;
    freqStop(u,:) = baseFrequency + timeZeroFreqOffset + chirp*scanTime/2;
end

end

function [opticalDist1,opticalDist2,opticalDist3] = SpecifyOpticalDistances()
    aodThickness = 8e-3;
    aodRefractiveIndexOrdOnAxis = 2.26;
    opticalDistOffset = aodThickness * (1 - 1/aodRefractiveIndexOrdOnAxis);

    distance1 = 50e-3;
    distance2 = 50e-3;
    distance3 = 50e-3;

    opticalDist1 = (distance1 + opticalDistOffset);
    opticalDist2 = (distance2 + opticalDistOffset);
    opticalDist3 = (distance3 + opticalDistOffset);
end

function y = ScaleAndRound(x,scaleFactor)
    y = round( x * scaleFactor );
end