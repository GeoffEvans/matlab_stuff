function [a, b, c, tRamp] = microscope_drivers_structural(...
    xyElements,...
    zElements,...
    acceptanceAngle,...
    dwellTimeMultiplier,...
    wavelength,...
    zoomInteger,...
    baseFrequency,...
    xCentre,...
    yCentre,...
    zCentre,...
    systemClockFreq)

V = 612.9;
aodAperture = 15e-3;
aodMode = -1;
dataTimeInterval = 50e-9;

normalisedCentres = CalculateNormalisedVoxelCentres(xyElements, zElements, zoomInteger);
[constantComponent, linearChirp, quadraticChirp] = CalculateFrequencyComponents();
[A, B, C] = ComputeReturnsForLabview();

    function normalisedCentres = CalculateNormalisedVoxelCentres(xyElements, zElements, zoomInteger)
        subdividedVox = xyElements * zoomInteger;
        [xGrid,yGrid,zGrid] = meshgrid(1:xyElements,1:xyElements,1:zElements);
        
        Mxgp = floor( Xp * subdividedVox / 2) + ceil(subdividedVox / 2); % zoomed voxel identifier
        Mygp = floor( Yp * subdividedVox / 2) + ceil(subdividedVox / 2);
        Mzgp = floor( Zp * subdividedVox / 2) + ceil(subdividedVox / 2);
        
        Mxp = Mxgp + xGrid - ceil(xyElements/2);
        Myp = Mygp + yGrid - ceil(xyElements/2);
        Mzp = Mzgp + zGrid - ceil(zElements/2);
        
        Xnp = (2*Mxp - subdividedVox + 1)/(subdividedVox - 1);
        Ynp = (2*Myp - subdividedVox + 1)/(subdividedVox - 1);
        Znp = (2*Mzp - subdividedVox + 1)/(subdividedVox - 1);
    end

    function [constantComponent, linearChirp, quadraticChirp] = CalculateFrequencyComponents()
        zCentre = zCentre + (zCentre == 0) * 1e-6;   % to prevent divide by zero
        % handle zero as a special case in if else
        
        opticalDist = OpticalDistance();
        
        aodFillTime = dataTimeInterval * ceil( aodAperture / V / dataTimeInterval ); % rounded up to a multiple
        halfFreqShiftMax = acceptanceAngle * V / wavelength;
        chirpMax = 2 * halfFreqShiftMax / aodFillTime;
        opticalDistMinimum4 = V^2 / (2 * wavelength * chirpMax);
        opticalDist(4) = opticalDistMinimum4 / zCentre; % ???
        
        chirp(1) = aodMode * (V^2/wavelength) / ( 2 * (opticalDist(1)+opticalDist(2)+opticalDist(3)+opticalDist(4)) );
        chirp(2) = aodMode * (V^2/wavelength) / ( 2 * (opticalDist(2)+opticalDist(3)+opticalDist(4)) );
        chirp(3) = aodMode * (V^2/wavelength) / ( 2 * (opticalDist(3)+opticalDist(4)) );
        chirp(4) = aodMode * (V^2/wavelength) / ( 2 * opticalDist(4) );
        
        timeZeroFreqOffset(1)=-aodMode*2*(opticalDist(3)+opticalDist(4))*opticalDist(4)/((opticalDist(4)+opticalDist(3))*(opticalDist(1)+opticalDist(2)+2*(opticalDist(3)+opticalDist(4))))*(Xnp(u)+skewZX.*Znp(u))*halfFreqShiftMax;%minus
        timeZeroFreqOffset(2)=aodMode*2*opticalDist(4)/(2*opticalDist(4)+opticalDist(3)+opticalDist(2))*(Ynp(u)+skewZY.*Znp(u))*halfFreqShiftMax;% AOD Y1 correct I think 150408
        timeZeroFreqOffset(3)=aodMode* opticalDist(4)/(opticalDist(4)+opticalDist(3))*(Xnp(u)+skewZX.*Znp(u))*halfFreqShiftMax;% AOD X2 correct I think 150408%minus
        timeZeroFreqOffset(4)=-aodMode* (Ynp(u)+skewZY.*Znp(u))*halfFreqShiftMax;% AOD Y2 correct I think 150408
        
        centreFreqOffset = timeZeroFreqOffset - chirp * aodFillTime/2;
    end

    function [A, B, C] = ComputeReturnsForLabview()
        A(u,:) = ScaleAndRound(baseFrequency + centreFreqOffset(1:4), 2^32 / systemClockFreq);
        B(u,:) = ScaleAndRound(chirp, 2^32 / systemClockFreq / 8); %CHANGED to scale B as for scanning mode
        C(u,:) = zero(numOfImageCentresPerCycle,4);
        
        dwellTime = dwellTimeMultiplier * dataTimeInterval;
        scanTime = aodFillTime + dwellTime; % the scan time for each point
        ticksPerRamp(u,:)=  repmat( ScaleAndRound(scanTime, systemClockFreq), 1, 4);
        freqStart(u,:) = baseFrequency + timeZeroFreqOffset - chirp*scanTime/2;
        freqStop(u,:) = baseFrequency + timeZeroFreqOffset + chirp*scanTime/2;
    end

    function opticalDist = OpticalDistance()
        aodThickness = 8e-3;
        aodRefractiveIndexOrdOnAxis = 2.26;
        opticalDistOffset = aodThickness * (1 - 1/aodRefractiveIndexOrdOnAxis);
        
        distance1 = 50e-3;
        distance2 = 50e-3;
        distance3 = 50e-3;
        
        opticalDist = zeros(4,1);
        opticalDist(1) = (distance1 + opticalDistOffset);
        opticalDist(2) = (distance2 + opticalDistOffset);
        opticalDist(3) = (distance3 + opticalDistOffset);
    end

    function y = ScaleAndRound(x,scaleFactor)
        y = round( x * scaleFactor );
    end
end

