function [baseFreqCompensated, linearChirp, rampTime, isMiniscan] = microscope_driver(...
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
    opWavelenVac)

isStructural = strcmp(imagingMode, 'structural');
isRaster = strcmp(imagingMode, 'raster');
isPointing = strcmp(imagingMode, 'pointing');
isMiniscan = strcmp(imagingMode, 'miniscan');

aodXyCentres = -1e-3 * [ [0;0], [2.3;0], [4.6;2.3], [4.6;4.6] ]; % negative sign is related to coordinate system and these are XY centers of the AODs to approx match the deflection shift.
referenceShift = [0;0;0];
commonModeOffsets = [4.65 5.45]; % commonMode is for eliminating scan speed displacement
differentialOffsets = [-0.85 3.6]; % differentialMde is for compensating XZ shifts

transducerOffsets = -1e-3 * [commonModeOffsets(1) - differentialOffsets(1) ...
                            commonModeOffsets(2) - differentialOffsets(2) ...
                            commonModeOffsets(1) + differentialOffsets(1) ... 
                            commonModeOffsets(2) + differentialOffsets(2)];

% for calculating the cartesian coordinates of the imaging region defined by the scan angle
referencePoint = ComputeReferencePoint(optimalFrequency, pairDeflectionRatio, aodMode, referenceShift, opWavelenVac); 

[xyDeflectionMm,focalLength,xyScanSpeed,rampTime] = ConvertNormalisedToCartesian(xyNumOfElems, zoomFactor, acceptanceAngle, imageCentreNormalised, scanDisplacementNormalised, referencePoint);

driveParams = GenerateDriveParams(xyDeflectionMm, pairDeflectionRatio, xyScanSpeed, optimalFrequency, focalLength, opWavelenVac);
[aodDirectionVectors, generalCentres, ~, aolDrives] = calculate_aol_drive(4, driveParams, aodMode);
[baseFreqCompensated, linearChirp] = CompensateFreqForTransducerLocation(aodXyCentres, transducerOffsets, generalCentres, aodDirectionVectors, aolDrives);

    function referencePoint = ComputeReferencePoint(optimalFrequency,pairDeflectionRatio,aodMode,referenceShift,opWavelenVac)
        [~, zeroDeflectionCentres, ~, ~] = calculate_aol_drive(4, aol_drive_params(inf, optimalFrequency, [0;0], pairDeflectionRatio, [0;0], opWavelenVac), aodMode);
        referencePoint = zeroDeflectionCentres{end} + referenceShift;
    end

    function z = RemoveNormalisationZ(zNormalised)
        z = aodAperture ./ (4 * acceptanceAngle * 1e-3 * zNormalised); % in m IF acceptanceAngle in mrad
    end

    function [xyDeflectionMm,zImageCentre,xyScanSpeed,rampTime] = ConvertNormalisedToCartesian(xyNumOfElems, zoomFactor, acceptanceAngle, imageCentreNorm, scanDispNorm, referencePoint)
        
        % Originally, half deflection went on each of the pair, so if max deflection on one is accAngle then total def is twice that, hence factor of 2 below
        if isRaster || isStructural
            zImageCentreNormalised = imageCentreNorm(3);
        else
            zImageCentreNormalised = imageCentreNorm(3,:); % one for each drive
        end
        zImageCentreNormalised = (zImageCentreNormalised == 0) * 1e-6 + zImageCentreNormalised; % avoid divide by 0
        zImageCentre = RemoveNormalisationZ(zImageCentreNormalised);
        
        xyExtremeRelToBaseRayMm = 2 * acceptanceAngle ./ zoomFactor .* zImageCentre; % in mm IF acceptanceAngle in mrad
        
        xyDeflectionMm = CalcXyDeflectionsForImagingMode(xyExtremeRelToBaseRayMm, xyNumOfElems, imageCentreNorm, referencePoint);
        
        [scanDispNorm, xyDeflectionMm, rampTime] = AlterDwellAndScanForImagingMode(scanDispNorm, xyDeflectionMm, referencePoint, xyNumOfElems, dwellTime);
        
        xScanSpeed = scanDispNorm(1,:) .* xyExtremeRelToBaseRayMm * 1e-3 ./ rampTime;
        yScanSpeed = scanDispNorm(2,:) .* xyExtremeRelToBaseRayMm * 1e-3 ./ rampTime;
        xyScanSpeed = [xScanSpeed; yScanSpeed];
        
        function xyDeflectionMm = CalcXyDeflectionsForImagingMode(xyExtremeRelToBaseRayMm, xyNumOfElems, imageCentreNormalised, referencePoint)

            xyImageCentreOffsetMm = [imageCentreNormalised(1,:) .* xyExtremeRelToBaseRayMm; imageCentreNormalised(2,:) .* xyExtremeRelToBaseRayMm];
            
            if isRaster
                xyImageCentreOffsetMm(1) = xyImageCentreOffsetMm(1) .* zoomFactor;
            end
            if isRaster || isStructural    
                
                xyRowRelToBaseRayMm = linspace( -xyExtremeRelToBaseRayMm, xyExtremeRelToBaseRayMm, xyNumOfElems);
                if xyNumOfElems == 1
                    xyRowRelToBaseRayMm = 0;
                end
                
                xMm = referencePoint(1)*1e3 + xyRowRelToBaseRayMm;% + xyImageCentreOffsetMm(1)*isStructural; % for now let's force structural and raster to focus at the same location...
                yMm = referencePoint(2)*1e3 + xyRowRelToBaseRayMm;% + xyImageCentreOffsetMm(2)*isStructural;
                
            else % must be in pointing or miniscan
                xMm = referencePoint(1)*1e3 + xyImageCentreOffsetMm(1,:);
                yMm = referencePoint(2)*1e3 + xyImageCentreOffsetMm(2,:);
            end
               
            xyDeflectionMm = [xMm;yMm];
        end
        
        function [scanDisplacementNormalised, xyDeflectionMm, rampTime] = AlterDwellAndScanForImagingMode(scanDisplacementNormalised, xyDeflectionMm, referencePoint, xyNumOfElems, dwellTime)
            if isRaster
                rampTime = dwellTime * xyNumOfElems;
                scanDisplacementNormalised = [2;0;0]; % from -1 to 1
                xyDeflectionMm(1,:) = 0*xyDeflectionMm(1,:) + referencePoint(1) * 1e3; % scans over x so centre at t=0
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

    function driveParams = GenerateDriveParams(xyDeflectionMm, pairDeflectionRatio, xyScanSpeed, optimalFrequency, focalLength, opWavelenVac)
        if isMiniscan || isPointing
            pdrs = stretch(pairDeflectionRatio, size(xyDeflectionMm,2));
            driveParams = aol_drive_params(focalLength, optimalFrequency, xyDeflectionMm, pdrs, xyScanSpeed, opWavelenVac);
        else
            driveParams = aol_drive_params.MakeDriveParams(xyDeflectionMm,pairDeflectionRatio,xyScanSpeed,optimalFrequency,focalLength,opWavelenVac);
        end
    end

    function [baseFreqCompensated, linearChirp] = CompensateFreqForTransducerLocation(aodXyCentres, transducerOffsets, baseRayCentres, aodDirectionVectors, aolDrives)
        baseFreq = [aolDrives.baseFreq];
        linearChirp = [aolDrives.chirp];
        
        timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, transducerOffsets, baseRayCentres, aodDirectionVectors);        

        baseFreqCompensated = baseFreq + linearChirp .* timeFromTransducerToBaseRay;
        
        function timeFromTransducerToBaseRay = CalculateTimeFromTransducerToBaseRay(aodXyCentres, transducerOffsets, baseRayCentres, aodDirectionVectors)
            numOfAods = numel(aodDirectionVectors);
            timeFromTransducerToBaseRay = zeros(numOfAods,size(baseRayCentres{1},2));
            for nthAod = 1:numOfAods
                xyDisplacementVector = Difference(baseRayCentres{nthAod}, aodXyCentres(:,nthAod));
                distanceFromTransducerToBaseRay = aodAperture/2 + transducerOffsets(nthAod) + XyDotProd(aodDirectionVectors{nthAod}, xyDisplacementVector);
                timeFromTransducerToBaseRay(nthAod,:) = distanceFromTransducerToBaseRay / V;
            end
            
            function dp = XyDotProd(u,v)
                dp = u(1)*v(1,:) + u(2)*v(2,:);
            end
            function diff = Difference(u,v)
                diff(1,:) = u(1,:) - v(1);
                diff(2,:) = u(2,:) - v(2);
            end
        end
    end
end

