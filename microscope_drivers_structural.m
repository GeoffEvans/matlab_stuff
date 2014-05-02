function [a, b, c, ticksPerRamp] = microscope_drivers_structural(...
    xyNumOfElems,...
    acceptanceAngle,... % mrad
    dwellTimeMultiplier,...
    wavelength,... % currently set in aod3d.m, can't imagine this ever changes
    zoom,...
    optimalFrequency,...
    xImageCentreNormalised,...
    yImageCentreNormalised,...
    zImageCentreNormalised,...
    systemClockFreq,...
    pairDeflectionRatio)

V = 612.8834; % As computed by teo2.find_v_ac_min(pi/2,pi/4)
aodAperture = 15e-3;
aodMode = -1;
dataTimeInterval = 50e-9;
xyImageCentreNormalised = [xImageCentreNormalised; yImageCentreNormalised];
aodXyCentres = [ [0;0], [0;0], [0;0], [0;0] ];

driveParams = aol_drive_params(1, optimalFrequency, [0;0], pairDeflectionRatio, 0);
[~, baseRayCentres, ~, ~] = calculate_aol_drive(4, driveParams);
[xyDeflectionMm,focalLength] = RemoveNormalisation(xyNumOfElems, zoom, acceptanceAngle, zImageCentreNormalised, xyImageCentreNormalised);
driveParams = aol_drive_params(focalLength, optimalFrequency, xyDeflectionMm, pairDeflectionRatio, 0);
[aodDirectionVectors, ~, ~, aolDrives] = calculate_aol_drive(4, driveParams);
aolDrives = CompensateFreqForTransducerLocation(aodXyCentres, baseRayCentres, aodDirectionVectors, aolDrives);
[a, b, c] = ComputeReturnsForLabview(aolDrives);

    function [xyDeflectionMm,zImageCentre] = RemoveNormalisation(xyNumOfElems, zoom, acceptanceAngle, zImageCentreNormalised, xyImageCentreNormalised)
        % Originally, half deflection went on each of the pair, so if max deflection on one is accAngle then total def is twice that, hence factor of 2 below
        zImageCentreNormalised = (zImageCentreNormalised == 0) * 1e-6 + zImageCentreNormalised; % avoid divide by 0
        zImageCentre = aodAperture / (4 * acceptanceAngle * zImageCentreNormalised);
        xyExtremeRelToBaseRay = 2 * acceptanceAngle / zoom * zImageCentre; % in mm assuming acceptanceAngle in mrad
        xRowRelToBaseRay = linspace( -xyExtremeRelToBaseRay, xyExtremeRelToBaseRay, xyNumOfElems); 
        [xGridRelToBaseRay,yGridRelToBaseRay] = meshgrid(xRowRelToBaseRay);

        xyImageCentre = xyImageCentreNormalised * xyExtremeRelToBaseRay;
        xGrid = xGridRelToBaseRay + baseRayCentres(1,4) + xyImageCentre(1);
        yGrid = yGridRelToBaseRay + baseRayCentres(2,4) + xyImageCentre(2);
        xyDeflectionMm = [xGrid(:); yGrid(:)]; 
    end

    function aolDrives = CompensateFreqForTransducerLocation(aodXyCentres, baseRayCentres, aodDirectionVectors, aolDrives)
        timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, baseRayCentres, aodDirectionVectors);
        aolDrives.baseFreq = aolDrives.baseFreq + aolDrives.chirp .* timeFromTransducerToBaseRay;

        function timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, baseRayCentres, aodDirectionVectors)
            distanceFromTransducerToBaseRay = aodAperture/2 + dot(aodDirectionVectors, baseRayCentres - aodXyCentres);
            timeFromTransducerToBaseRay = distanceFromTransducerToBaseRay / V;
        end
    end

    function [a, b, c] = ComputeReturnsForLabview(aolDrives)

        a = ScaleAndRound(aolDrives.baseFreq, aodMode * 2^32 / systemClockFreq);
        b = ScaleAndRound(aolDrives.chirp, aodMode * 2^32 / systemClockFreq / 8); % scale B down here and expand later
        c = b * 0;
        
        dwellTime = dwellTimeMultiplier * dataTimeInterval;
        scanTime = aodFillTime + dwellTime; % the scan time for each point
        ticksPerRamp =  repmat( ScaleAndRound(scanTime, systemClockFreq), 1, 4);
    end

    function y = ScaleAndRound(x,scaleFactor)
        y = round( x * scaleFactor );
    end
end

