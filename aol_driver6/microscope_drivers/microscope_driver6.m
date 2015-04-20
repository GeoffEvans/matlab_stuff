function [baseFreqCompensated, linearChirp, rampTime, isMiniscan] = microscope_driver6(...
    imagingMode,...    
    aodMode,...
    xyNumOfElems,...
    acceptanceAngle,... % mrad
    dwellTime,...
    zoomFactor,...
    optimalFrequency,...
    imageCentreNormalised,...
    scanDisplacementNormalised,...
    aodAperture,...
    V,...
    opWavelenVac)

if aodMode ~= -1
    error('aolMode must be -1')
end 
if any(scanDisplacementNormalised(2:3,:) ~= 0)
    error('only scan in x')
end 

isStructural = strcmp(imagingMode, 'structural');
isRaster = strcmp(imagingMode, 'raster');
isPointing = strcmp(imagingMode, 'pointing');
isMiniscan = strcmp(imagingMode, 'miniscan');

aodXyCentres = -1e-3 * [ [0;0], [0;0], [0;0], [0;0], [0;0], [0;0] ]; % negative sign is because first AOD deflects in -x direction
referenceShiftMm = [0;0;0];

transducerOffsets = -1e-3 * [0,0,0,0,0,0];

[xyDeflectionMm,focalLength,xyScanSpeed,rampTime] = ConvertNormalisedToCartesian(xyNumOfElems, zoomFactor, acceptanceAngle, imageCentreNormalised, scanDisplacementNormalised, referenceShiftMm);
driveParams = GenerateDriveParams(xyDeflectionMm, xyScanSpeed, optimalFrequency, focalLength, opWavelenVac);
[baseRayCentres, aodDirectionVectors, ~, aolDrives] = calculate_aol_drive6(driveParams);
[baseFreqCompensated, linearChirp] = CompensateFreqForTransducerLocation(aodXyCentres, transducerOffsets, baseRayCentres, aodDirectionVectors, aolDrives, V);

    function z = RemoveNormalisationZ(zNormalised, acceptanceAngle)
        z = aodAperture ./ (4 * acceptanceAngle * 1e-3 * zNormalised); % in m IF acceptanceAngle in mrad
    end

    function [xyDeflectionMm,zImageCentre,xyScanSpeed,rampTime] = ConvertNormalisedToCartesian(xyNumOfElems, zoomFactor, acceptanceAngle, imageCentreNorm, scanDispNorm, referenceShiftMm)
        if isRaster || isStructural
            zImageCentreNormalised = imageCentreNorm(3);
        else
            zImageCentreNormalised = imageCentreNorm(3,:); % one for each drive
        end
        zImageCentreNormalised = (zImageCentreNormalised == 0) * 1e-6 + zImageCentreNormalised; % avoid divide by 0
        
        zImageCentre = RemoveNormalisationZ(zImageCentreNormalised, acceptanceAngle) + referenceShiftMm(3)*1e-3;
        xyExtremeRelToBaseRayMm = 2 * acceptanceAngle ./ zoomFactor .* zImageCentre; % in mm IF acceptanceAngle in mrad % Originally, half deflection went on each of the pair, so if max deflection on one is accAngle then total def is twice that, hence factor of 2 
        xyDeflectionMm = CalcXyDeflectionsForImagingMode(xyExtremeRelToBaseRayMm, xyNumOfElems, imageCentreNorm, referenceShiftMm, zoomFactor);
        
        [scanDispNorm, xyDeflectionMm, rampTime] = AlterDwellAndScanForImagingMode(scanDispNorm, xyDeflectionMm, xyNumOfElems, dwellTime);
        
        xScanSpeed = scanDispNorm(1,:) .* xyExtremeRelToBaseRayMm * 1e-3 ./ rampTime;
        yScanSpeed = scanDispNorm(2,:) .* xyExtremeRelToBaseRayMm * 1e-3 ./ rampTime;
        xyScanSpeed = [xScanSpeed; yScanSpeed];
        
        function xyDeflectionMm = CalcXyDeflectionsForImagingMode(xyExtremeRelToBaseRayMm, xyNumOfElems, imageCentreNormalised, referenceShiftMm, zoomFactor)

            xyImageCentreOffsetMm = [imageCentreNormalised(1,:) .* xyExtremeRelToBaseRayMm; imageCentreNormalised(2,:) .* xyExtremeRelToBaseRayMm];
            
            if isRaster
                xyImageCentreOffsetMm(1) = xyImageCentreOffsetMm(1) .* zoomFactor;
            end
            if isRaster || isStructural    
                
                xyRowRelToBaseRayMm = linspace( -xyExtremeRelToBaseRayMm, xyExtremeRelToBaseRayMm, xyNumOfElems);
                if xyNumOfElems == 1
                    xyRowRelToBaseRayMm = 0;
                end
                
                xMm = xyRowRelToBaseRayMm;% + xyImageCentreOffsetMm(1)*isStructural; % for now let's force structural and raster to focus at the same location...
                yMm = xyRowRelToBaseRayMm;% + xyImageCentreOffsetMm(2)*isStructural;
                
            else % must be in pointing or miniscan
                xMm = xyImageCentreOffsetMm(1,:);
                yMm = xyImageCentreOffsetMm(2,:);
            end
               
            xyDeflectionMm = [xMm + referenceShiftMm(1); yMm + referenceShiftMm(2)];
        end
        
        function [scanDisplacementNormalised, xyDeflectionMm, rampTime] = AlterDwellAndScanForImagingMode(scanDisplacementNormalised, xyDeflectionMm, xyNumOfElems, dwellTime)
            if isRaster
                rampTime = dwellTime * xyNumOfElems;
                scanDisplacementNormalised = [2;0;0]; % from -1 to 1
                xyDeflectionMm(1,:) = 0*xyDeflectionMm(1,:); % scans over x so centre at t=0
            elseif isStructural
                rampTime = dwellTime;
                scanDisplacementNormalised = [0;0;0];
            elseif isPointing
                rampTime = dwellTime;
                scanDisplacementNormalised = scanDisplacementNormalised .* 0;
            elseif isMiniscan
                rampTime = dwellTime * xyNumOfElems * mag(scanDisplacementNormalised)/2;
                % miniscan keeps defaults
            else
                error('imaging mode invalid')
            end
        end
    end

    function driveParams = GenerateDriveParams(xyDeflectionMm, xyScanSpeed, optimalFrequency, focalLength, opWavelenVac)
        if isMiniscan || isPointing
            driveParams = aol_drive_params6(focalLength, optimalFrequency, xyDeflectionMm, xyScanSpeed, opWavelenVac);
        else
            driveParams = aol_drive_params6.MakeDriveParams(xyDeflectionMm,xyScanSpeed,optimalFrequency,focalLength,opWavelenVac);
        end
    end

    function [baseFreqCompensated, linearChirp] = CompensateFreqForTransducerLocation(aodXyCentres, transducerOffsets, baseRayCentres, aodDirectionVectors, aolDrives, V)
        baseFreq = [aolDrives.baseFreq];
        linearChirp = [aolDrives.chirp];
        
        timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, transducerOffsets, baseRayCentres, aodDirectionVectors, V);        

        baseFreqCompensated = baseFreq + linearChirp .* transpose(timeFromTransducerToBaseRay);
        
        function timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, transducerOffsets, baseRayCentres, aodDirectionVectors, V)
            numOfAods = size(aodDirectionVectors,2);
            timeFromTransducerToBaseRay = zeros(numOfAods,size(baseRayCentres,3));
            for nthAod = 1:numOfAods
                xyDisplacementVector = squeeze([baseRayCentres(1,nthAod,:) - aodXyCentres(1,nthAod); baseRayCentres(2,nthAod,:) - aodXyCentres(2,nthAod)]);
                distanceFromTransducerToBaseRay = aodAperture/2 + transducerOffsets(nthAod) + XyDotProd(aodDirectionVectors(:,nthAod), xyDisplacementVector);
                timeFromTransducerToBaseRay(nthAod,:) = distanceFromTransducerToBaseRay / V;
            end
            
            function dp = XyDotProd(u,v)
                dp = u(1)*v(1,:) + u(2)*v(2,:);
            end
        end
    end
end

