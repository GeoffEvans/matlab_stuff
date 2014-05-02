function [a, b, c, ticksPerRamp] = microscope_drivers_structural(...
    xyNumOfElems,...
    acceptanceAngle,... % mrad
    dwellTimeMultiplier,...
    wavelength,... % currently set in aod3d.m, can't imagine this ever changes
    zoom,...
    optimalFrequency,...
    xImageCentre,...
    yImageCentre,...
    zImageCentre,...
    systemClockFreq,...
    pairDeflectionRatio)

V = 612.8834; % As computed by teo2.find_v_ac_min(pi/2,pi/4)
aodAperture = 15e-3;
aodMode = -1;
dataTimeInterval = 50e-9;
xyImageCentre = [xImageCentre; yImageCentre];
aodXyCentres = [ [0;0], [0;0], [0;0], [0;0] ];

xyDeflectionMm = SpecifyXyGrid(xyNumOfElems, zoom, acceptanceAngle, zImageCentre, xyImageCentre);
driveParams = aol_drive_params(focalLength, optimalFrequency, xyDeflectionMm, pairDeflectionRatio, 0);
[aodDirectionVectors, baseRayCentres, ~, aolDrives] = calculate_aol_drive(4, driveParams);
aolDrives = CompensateFreqForTransducerLocation(aodXyCentres, baseRayCentres, aodDirectionVectors, aolDrives);
[a, b, c] = ComputeReturnsForLabview(aolDrives);

    function xyDeflectionMm = SpecifyXyGrid(xyNumOfElems, zoom, acceptanceAngle, zImageCentre, xyImageCentre)
        % Originally, half deflection went on each of the pair, so if max deflection on one is accAngle 
        % then total def is twice that, hence factor of 2 below
        xyExtremeRelToBaseRay = 2 * acceptanceAngle / zoom * zImageCentre + xyImageCentre; % in mm assuming acceptanceAngle in mrad
        xRowRelToBaseRay = linspace( -xyExtremeRelToBaseRay, xyExtremeRelToBaseRay, xyNumOfElems); 
        [xGridRelToBaseRay,yGridRelToBaseRay] = meshgrid(xRowRelToBaseRay);
        xGrid = xGridRelToBaseRay + baseRayCentres(1,4);
        yGrid = yGridRelToBaseRay + baseRayCentres(2,4);
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

