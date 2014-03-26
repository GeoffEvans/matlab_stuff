function [aodDirectionVectors, aodCentres, zFocusPredicted, aolDrives] = calculate_aol_drive(numOfAods, driveParams)
% Returns cell aodDirectionVectors; 3,numOfAods aodCentres; numOfAods,numOfDrives chirp; numOfAods,numOfDrives baseFreq}

optimalBaseFreq = driveParams.optimalBaseFreq;
xyDeflectionMm = driveParams.xyDeflectionMm;
pairDeflectionRatio = driveParams.pairDeflectionRatio;
scanSpeed = driveParams.scanSpeed;
focalLength = driveParams.focalLength;

widthOfAod = 5e-3;
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
[aodCentres,zFocusPredicted] = CalculateAodCentres(numOfAods,optimalBaseFreq,aodSpacing,aodDirectionVectors); % use optimalBaseFreq because this needs to be fixed over pointing fov.
aolDrives = aol_drive_freqs(baseFreq,chirp);


    function [centres,zFocusPredicted] = CalculateAodCentres(numOfAods,optimalBaseFreq,aodSpacing,aodDirectionVectors)
        centres = zeros(3,numOfAods);% calc centre of AODs and the focal point using small angle approx.
        
        for aodNumber = 2:numOfAods % start from 2 because offset on first AOD is zero
            centres(1:2,aodNumber) = XyDeflectionThroughPriorAods(aodNumber,optimalBaseFreq,aodSpacing,aodDirectionVectors);
        end
        
        zPlanesAod = 1 + [0, cumsum(aodSpacing)];
        centres(3,:) = zPlanesAod(1:end-1);
        zFocusPredicted = zPlanesAod(end);
        
        function xyDeflectionCalc = XyDeflectionThroughPriorAods(aodNumber,optimalBaseFreq,aodSpacing,aodDirectionVectors)
            xyDeflectionCalc = 0; % initially zero deflection
            for mthDeflectingAod = 1:aodNumber-1
                offsetDueToMthAod = sum(aodSpacing(mthDeflectingAod:(aodNumber-1))).*aodDirectionVectors{mthDeflectingAod}.*optimalBaseFreq; 
                xyDeflectionCalc = xyDeflectionCalc + offsetDueToMthAod;
            end
            xyDeflectionCalc = xyDeflectionCalc * lambda/V;
        end
    end

    function [aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = Aod4(focalLength, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed)
        %aodDirectionVectors = {[0;1], [1;0], -[0;1], -[1;0]};
        aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
        aodSpacing = [5e-2, 5e-2, 5e-2];
        aodSpacing = [aodSpacing, sum(aodSpacing)+focalLength];
        
        effctvSpcng = aodSpacing - correctionDistance;
        
        A = scanSpeed / V;
        chirpFactor = [(1 + A)./(effctvSpcng(1) + effctvSpcng(2) + 2*effctvSpcng(3) + 2*effctvSpcng(4) + A*effctvSpcng(1) + A*effctvSpcng(2));...
            A*0 + 1./(effctvSpcng(2) + effctvSpcng(3) + 2*effctvSpcng(4));...
            (1 - A)./(2*(effctvSpcng(3) + effctvSpcng(4)));...
            A*0 + 1./(2*effctvSpcng(4))];
        
        f3diff = - V/lambda .* xyDeflectionMm(1,:) .* 1e-3 ./ (pairDeflectionRatio .* sum(effctvSpcng) + sum(effctvSpcng(3:4)));
        f4diff = - V/lambda .* xyDeflectionMm(2,:) .* 1e-3 ./ (pairDeflectionRatio .* sum(effctvSpcng(2:4)) + effctvSpcng(4));
        baseFreq = [ sum(effctvSpcng(3:4))/sum(effctvSpcng).*optimalBaseFreq - pairDeflectionRatio.*f3diff;...
            effctvSpcng(4)/sum(effctvSpcng(2:4)).*optimalBaseFreq - pairDeflectionRatio.*f4diff; optimalBaseFreq+f3diff; optimalBaseFreq+f4diff];
    end
end