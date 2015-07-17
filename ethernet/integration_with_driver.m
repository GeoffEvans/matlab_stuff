function [a, b, ticksPerRamp] = integration_with_driver()
    opWavelenVac = 850 * 1e-9;
    optimalFrequency = [1,1,1,1] * 35 * 1e6;
    acceptanceAngle = 3.5; % mrad
    zoomFactor = 1;
    xyNumOfElems = 100;
    dwellTime = 1e-6 * 1;
    imageCentreNormalised = [0;0;0]; 
    scanDisplacementNormalised = [0;0;0];
    pairDeflectionRatio = -1;
    aodAperture = 15e-3;
    V = 612.8;
    aodMode = -1;
    imagingMode = 'structural';
    systemClockFreq = 200e6;
    dataTimeInterval = 1/systemClockFreq;
    aodAperture = 16e-3;
    
    [baseFreq, linearChirp, rampTime, isMiniscan] = microscope_driver(...
        imagingMode,...    
        aodMode,...
        xyNumOfElems,...
        acceptanceAngle,... % mrad
        dwellTime,...
        zoomFactor,...
        optimalFrequency,...
        imageCentreNormalised,...
        scanDisplacementNormalised,...
        pairDeflectionRatio,...
        aodAperture,...
        V,...
        opWavelenVac);

    a = SwapScaleAndRound(baseFreq', 2^32 / systemClockFreq);
    b = SwapScaleAndRound(linearChirp', 2^32 / systemClockFreq^2 * 2 / 2^3); % scale B down here and expand later

    aodFillTime = dataTimeInterval * ceil(aodAperture / V / dataTimeInterval);
    scanTime = aodFillTime + rampTime; % the scan time for each point

    if isMiniscan
        ticksPerRamp =  SwapScaleAndRound( repmat(scanTime(:), 1, 4), systemClockFreq);
    else
        ticksPerRamp =  SwapScaleAndRound( repmat(scanTime, 200, 4), systemClockFreq);
    end

    function y = SwapScaleAndRound(x,scaleFactor)
        y = round( x * scaleFactor );
        y(:,[3 2]) = y(:,[2 3]);
    end
end