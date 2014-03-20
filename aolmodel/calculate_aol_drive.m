function [aodDirectionVectors, aodCentres, zFocusPredicted, aolDrives] = calculate_aol_drive(numOfAods, optimalBaseFreq, xyDeflectionMm, pairDeflectionRatio, scanSpeed)
% Returns 2xnumOfAods aodDirectionVectors, 3xnumOfAods aodCentres, {numOfDrivesxnumOfAods chirp, numOfDrivesxnumOfAods baseFreq}
if length(optimalBaseFreq) ~= 1
    error('optimal base freq must be singleton')
end

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

[xyDeflectionMm, pairDeflectionRatio, scanSpeed] = GetAllCombinations(xyDeflectionMm, pairDeflectionRatio, scanSpeed);
[aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = AodDriveFunction(xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed);

chirp = V*V/lambda * chirpFactor;
[aodCentres,zFocusPredicted] = CalculateAodCentres(numOfAods,optimalBaseFreq,aodSpacing,aodDirectionVectors); % use optimalBaseFreq because this needs to be fixed over pointing fov.
aolDrives = aol_drive(baseFreq,chirp);


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

    function [aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = Aod4(xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed)
        %aodDirectionVectors = {[0;1], [1;0], -[0;1], -[1;0]};
        aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
        aodSpacing = [5e-2, 5e-2, 5e-2, 2];
        
        l1 = aodSpacing(1) - correctionDistance;
        l2 = aodSpacing(2) - correctionDistance;
        l3 = aodSpacing(3) - correctionDistance;
        l4 = aodSpacing(4) - correctionDistance;
        
        A = scanSpeed / V;
        chirpFactor = [(1 + A)/(l1 + l2 + 2*l3 + 2*l4 + A*l1 + A*l2)...
            1/(l2 + l3 + 2*l4)...
            (1 - A)/(2*(l3 + l4))...
            1/(2*l4)];
        
        f3diff = - V/lambda .* xyDeflectionMm(1,:) .* 1e-3 ./ (pairDeflectionRatio .* sum(aodSpacing) + sum(aodSpacing(3:4)));
        f4diff = - V/lambda .* xyDeflectionMm(2,:) .* 1e-3 ./ (pairDeflectionRatio .* sum(aodSpacing(2:4)) + aodSpacing(4));
        baseFreq = [ sum(aodSpacing(3:4))/sum(aodSpacing).*optimalBaseFreq - pairDeflectionRatio.*f3diff,...
            aodSpacing(4)/sum(aodSpacing(2:4)).*optimalBaseFreq - pairDeflectionRatio.*f4diff, optimalBaseFreq+f3diff, optimalBaseFreq+f4diff];
    end
end

function [twoByMany,b,c] = GetAllCombinations(twoByMany,b,c)
    [b,c] = meshgrid(b,c);
    b = b(:);
    c = c(:);
    twoByMany = stretch(twoByMany,length(b));
end