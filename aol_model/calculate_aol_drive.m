function [aodDirectionVectors, aodCentres, zFocusPredicted, aolDrives] = calculate_aol_drive(numOfAods, driveParams)
% Returns (cell()) aodDirectionVectors; (cell[numOfAods](3,numOfDrives)) aodCentres; (numOfAods,numOfDrives) chirp; numOfAods,numOfDrives baseFreq}

optimalBaseFreq = driveParams.optimalBaseFreq;
xyDeflectionMm = driveParams.xyDeflectionMm;
pairDeflectionRatio = driveParams.pairDeflectionRatio;
scanSpeed = driveParams.scanSpeed;
focalLength = driveParams.focalLength;

widthOfAod = 8e-3;
correctionDistance = widthOfAod * ( 1 - 1/2.26 );
V = teo2.find_v_ac_min(pi/2,pi/4);
lambda = aod3d.opWavelenVac;

AodDriveFunction = @Aod4;
switch(numOfAods)
   case 1 
      AodDriveFunction = @Aod1;
   case 2 
      AodDriveFunction = @Aod2;
end

[aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = AodDriveFunction(focalLength, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed);

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
            for mthDeflectingAod = 1:aodNumber-1
                spacings = aodSpacing(mthDeflectingAod:(aodNumber-1));
                offsetDueToMthAod = sum(spacings) .* aodDirectionVectors{mthDeflectingAod} * baseFreq(mthDeflectingAod,:); % note use of outer product here
                xyDeflectionCumulator = xyDeflectionCumulator + offsetDueToMthAod;
            end
            xyDeflectionCumulator = xyDeflectionCumulator * lambda/V;
        end
    end

    function [aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = Aod4(focalLength, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed)
        %aodDirectionVectors = {[0;1], [1;0], -[0;1], -[1;0]};
        aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
        aodSpacing = [5e-2, 5e-2, 5e-2, focalLength];
        
        effctvSpcng = aodSpacing - correctionDistance;
        
        A = scanSpeed / V;
        chirpFactor = [(1 + A)./(effctvSpcng(1) + effctvSpcng(2) + 2*effctvSpcng(3) + 2*effctvSpcng(4) + A*effctvSpcng(1) + A*effctvSpcng(2));...
            A*0 + 1./(effctvSpcng(2) + effctvSpcng(3) + 2*effctvSpcng(4));...
            (1 - A)./(2*(effctvSpcng(3) + effctvSpcng(4)));...
            A*0 + 1./(2*effctvSpcng(4))];
        
        % for constant components, choose r to represent the ratio of ANGULAR deflection on the first of the pair to the second of the pair
        % traditionally, this means r = 1, while all on the second would be r = 0
        dfx = (V/lambda .* xyDeflectionMm(1,:).*1e-3 - optimalBaseFreq .* sum(effctvSpcng(1:2))) ./ (pairDeflectionRatio .* sum(effctvSpcng) + sum(effctvSpcng(3:4)));
        dfy = (V/lambda .* xyDeflectionMm(2,:).*1e-3 - optimalBaseFreq .* sum(effctvSpcng(2:3))) ./ (pairDeflectionRatio .* sum(effctvSpcng(2:4)) + effctvSpcng(4));
        baseFreq = [ optimalBaseFreq + pairDeflectionRatio.*dfx; optimalBaseFreq + pairDeflectionRatio.*dfy; optimalBaseFreq - dfx; optimalBaseFreq - dfy ];
    end
end