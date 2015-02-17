function [aodDirectionVectors, aodCentres, zFocusPredicted, aolDrives] = calculate_aol_drive(numOfAods, driveParams, aodMode)
% Returns (cell()) aodDirectionVectors; (cell[numOfAods](3,numOfDrives)) aodCentres; (numOfAods,numOfDrives) chirp; numOfAods,numOfDrives baseFreq}

optimalBaseFreq = driveParams.optimalBaseFreq;
xyDeflectionMm = driveParams.xyDeflectionMm;
pairDeflectionRatio = driveParams.pairDeflectionRatio;
xyScanSpeed = driveParams.xyScanSpeed;
focalLengths = driveParams.focalLength;
lambda = driveParams.opWavelenVac;

widthOfAod = 8e-3;
correctionDistance = 0; %widthOfAod * ( 1 - 1/2.26 );
V = teo2.find_v_ac_min(pi/2,pi/4);

[aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = Aod4(focalLengths, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, xyScanSpeed);

chirp = V*V/lambda * chirpFactor;
[aodCentres,zFocusPredicted] = CalculateAodCentres(numOfAods,baseFreq,aodSpacing,aodDirectionVectors); % use optimalBaseFreq because this needs to be fixed over pointing fov.
aolDrives = aol_drive_freqs(baseFreq,chirp);

    function [aodCentresCell,zFocusPredicted] = CalculateAodCentres(numOfAods,baseFreq,aodSpacing,aodDirectionVectors)
        numOfDrives = size(baseFreq,2);
        
        zPlanesAod = 1 + [0, cumsum(aodSpacing)];
        zFocusPredicted = zPlanesAod(end);
        
        aodCentresCell = cell(1,numOfAods); % calc centre of AODs and the focal point using small angle approx.
        aodCentresCell{1} = [zeros(2,numOfDrives); ones(1,numOfDrives) * zPlanesAod(1)]; % first aod centre defines xy-origin
        
        for aodNumber = 2:numOfAods % start from 2 because offset on first AOD is zero
            xyCentres = XyDeflectionThroughPriorAods(aodNumber,numOfDrives,baseFreq,aodSpacing,aodDirectionVectors);
            zCentre = zPlanesAod(aodNumber);
            aodCentresCell{aodNumber} = [xyCentres; ones(1,numOfDrives) * zCentre]; 
        end
        
        function xyDeflectionCumulator = XyDeflectionThroughPriorAods(aodNumber,numOfDrives,baseFreq,aodSpacing,aodDirectionVectors)
            xyDeflectionCumulator = zeros(2,numOfDrives);
            for mthAod = 1:aodNumber-1
                effctvSpacesFromMthToAodNum = aodSpacing(mthAod:(aodNumber-1)) - correctionDistance;
                offsetDueToMthAod = sum(effctvSpacesFromMthToAodNum) .* aodDirectionVectors{mthAod} * baseFreq(mthAod,:); % note use of outer product here
                xyDeflectionCumulator = xyDeflectionCumulator + offsetDueToMthAod;
            end
            xyDeflectionCumulator = xyDeflectionCumulator * aodMode*lambda/V;
        end
    end

    function [aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = Aod4(focalLengths, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, xyScanSpeed)
        aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
        aodSpacing = [5e-2, 5e-2, 5e-2];
        
        effctvSpcng = aodSpacing - correctionDistance;
        effctvFcs = focalLengths - correctionDistance; % seperate to handle many focal lengths
        
        A = xyScanSpeed / V;
        chirpFactor = [(1 + A(1,:))./(effctvSpcng(1) + effctvSpcng(2) + 2*effctvSpcng(3) + 2*effctvFcs + A(1,:)*effctvSpcng(1) + A(1,:)*effctvSpcng(2));...
            (1 + A(2,:))./(effctvSpcng(2) + effctvSpcng(3) + 2*effctvFcs + A(2,:)*effctvSpcng(2) + A(2,:)*effctvSpcng(3));... 
            (1 - A(1,:))./(2 * (effctvSpcng(3) + effctvFcs) );...
            (1 - A(2,:))./(2 * effctvFcs)] / aodMode;
        
        % for constant components, choose r to represent the ratio of ANGULAR deflection on the first of the pair to the second of the pair
        % traditionally, this means r = 1, while all on the second would be r = 0
        dfx = ( V/lambda/aodMode .* xyDeflectionMm(1,:).*1e-3 - optimalBaseFreq(1) .* (sum(effctvSpcng) + effctvFcs) + optimalBaseFreq(3) .* (effctvSpcng(3) + effctvFcs) ) ...
                ./ ( pairDeflectionRatio .* (sum(effctvSpcng) + effctvFcs) + (effctvSpcng(3) + effctvFcs) );
        dfy = ( V/lambda/aodMode .* xyDeflectionMm(2,:).*1e-3 - optimalBaseFreq(2) .* (sum(effctvSpcng(2:3)) + effctvFcs) + optimalBaseFreq(4) .* effctvFcs ) ...
                ./ ( pairDeflectionRatio .* (sum(effctvSpcng(2:3)) + effctvFcs) + effctvFcs );
        dfx(focalLengths == inf) = 0; % for calculating coordinate reference points
        dfy(focalLengths == inf) = 0;
        baseFreq = [ optimalBaseFreq(1) + pairDeflectionRatio.*dfx; optimalBaseFreq(2) + pairDeflectionRatio.*dfy; optimalBaseFreq(3) - dfx; optimalBaseFreq(4) - dfy ];
    end
end