function [a, b, c, ticksPerRamp] = microscope_drivers_structural(...
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
    pairDeflectionRatio)

[V, aodAperture, aodMode, dataTimeInterval ] = DeclareConsts();
aodXyCentres = [ [0;0], [0;0], [0;0], [0;0] ];
xyImageCentreNormalised = [xImageCentreNormalised; yImageCentreNormalised];

[~, zeroDeflectionCentres, ~, ~] = calculate_aol_drive(4, aol_drive_params(inf, optimalFrequency, [0;0], pairDeflectionRatio, 0));

[xyDeflectionMm,focalLength] = RemoveNormalisation(xyNumOfElems, zoomFactor, acceptanceAngle, zImageCentreNormalised, xyImageCentreNormalised, zeroDeflectionCentres{end});

[aodDirectionVectors, generalCentres, ~, aolDrives] = calculate_aol_drive(4, aol_drive_params.MakeDriveParams(xyDeflectionMm,pairDeflectionRatio,0,optimalFrequency,focalLength));

[baseFreqCompensated,chirp] = CompensateFreqForTransducerLocation(aodXyCentres, zeroDeflectionCentres, aodDirectionVectors, aolDrives);

[a, b, c] = ComputeReturnsForLabview(baseFreqCompensated,chirp);

    function [xyDeflectionMm,zImageCentre] = RemoveNormalisation(xyNumOfElems, zoomFactor, acceptanceAngle, zImageCentreNormalised, xyImageCentreNormalised, xyzFinalAodRayCentres)
        
        % Originally, half deflection went on each of the pair, so if max deflection on one is accAngle then total def is twice that, hence factor of 2 below
        zImageCentreNormalised = (zImageCentreNormalised == 0) * 1e-6 + zImageCentreNormalised; % avoid divide by 0
        zImageCentre = aodAperture / (4 * acceptanceAngle * zImageCentreNormalised); % in m
        
        xyExtremeRelToBaseRay = 2 * acceptanceAngle / zoomFactor * zImageCentre; % in mm assuming acceptanceAngle in mrad
        
        xRowRelToBaseRay = linspace( -xyExtremeRelToBaseRay, xyExtremeRelToBaseRay, xyNumOfElems);
        [xGridRelToBaseRay,yGridRelToBaseRay] = meshgrid(xRowRelToBaseRay);
        
        xyImageCentreOffset = xyImageCentreNormalised * xyExtremeRelToBaseRay;
        
        xGrid = xGridRelToBaseRay + xyzFinalAodRayCentres(1) + xyImageCentreOffset(1);
        yGrid = yGridRelToBaseRay + xyzFinalAodRayCentres(2) + xyImageCentreOffset(2);
        xyDeflectionMm = [xGrid(:), yGrid(:)]';
    end

    function [baseFreqCompensated,chirp] = CompensateFreqForTransducerLocation(aodXyCentres, baseRayCentres, aodDirectionVectors, aolDrives)
        baseFreq = reshape([aolDrives.baseFreq],4,numel(aolDrives));
        chirp = reshape([aolDrives.chirp],4,numel(aolDrives));
        
        timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, cell2mat(baseRayCentres), cell2mat(aodDirectionVectors));        

        baseFreqCompensated = baseFreq + chirp .* repmat(timeFromTransducerToBaseRay(:),1,numel(aolDrives)); % QQ TODO amend to use baseRayCentre for each drive
        
        function timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, baseRayCentres, aodDirectionVectors)
            distanceFromTransducerToBaseRay = aodAperture/2 + dot(aodDirectionVectors, baseRayCentres(1:2,:) - aodXyCentres);
            timeFromTransducerToBaseRay = distanceFromTransducerToBaseRay / V;
        end
    end

    function [a, b, c] = ComputeReturnsForLabview(baseFreq,chirp)
        
        a = ScaleAndRound(baseFreq, aodMode * 2^32 / systemClockFreq);
        b = ScaleAndRound(chirp, aodMode * 2^32 / systemClockFreq / 8); % scale B down here and expand later
        c = b * 0;
        
        dwellTime = dwellTimeMultiplier * dataTimeInterval;
        aodFillTime = aodAperture / V;
        scanTime = aodFillTime + dwellTime; % the scan time for each point
        ticksPerRamp =  repmat( ScaleAndRound(scanTime, systemClockFreq), 1, 4);
    end


end

function y = ScaleAndRound(x,scaleFactor)
y = round( x * scaleFactor );
end

function [V, aodAperture, aodMode, dataTimeInterval ] = DeclareConsts()
V = 612.8834; % As computed by teo2.find_v_ac_min(pi/2,pi/4)
aodAperture = 15e-3;
aodMode = -1;
dataTimeInterval = 50e-9;
end