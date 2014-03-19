function [aodDirectionVectors, aodCentresAndFocus, aolDrives] = calculate_aol_drive(numOfAods, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed)
% Returns 2xnumOfAods aodDirectionVectors, 3xnumOfAods aodCentres, {numOfDrivesxnumOfAods chirp, numOfDrivesxnumOfAods baseFreq}

widthOfAod = 5e-3;
correctionDistance = widthOfAod * ( 1 - 1/2.26 );
V = teo2.find_v_ac_min(pi/2,pi/4);
lambda = aod3d.opWavelenVac;

switch(numOfAods)
   case 1 
      AodDriveFunction = @Aod1;
   case 2 
      AodDriveFunction = @Aod2;
   case 4 
      AodDriveFunction = @Aod4;
end

[aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = AodDriveFunction(xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed);

chirp = V*V/lambda * chirpFactor;
aodCentresAndFocus = CalculateAodCentres(numOfAods,optimalBaseFreq,aodSpacing,aodDirectionVectors); % use optimalBaseFreq because this needs to be fixed over pointing fov.
aolDrives = aol_drive(baseFreq,chirp);

    function centresAndFocus = CalculateAodCentres(numOfAods,optimalBaseFreq,aodSpacing,aodDirectionVectors)
        centresAndFocus = zeros(3,numOfAods+1);% calc centre of AODs and the focal point using small angle approx.
        
        for aodNumber = 2:numOfAods+1 % start from 2 because offset on first AOD is zero
            centresAndFocus(1:2,aodNumber) = XyDeflectionThroughPriorAods(aodNumber,optimalBaseFreq,aodSpacing,aodDirectionVectors);
        end
        
        zPlanesAod = 1 + [0, cumsum(aodSpacing)];
        centresAndFocus(3,:) = zPlanesAod;
        
        function xyDeflection = XyDeflectionThroughPriorAods(aodNumber,optimalBaseFreq,aodSpacing,aodDirectionVectors)
            xyDeflection = 0; % initially zero deflection
            for mthDeflectingAod = 1:aodNumber-1
                offsetDueToMthAod = sum(aodSpacing(mthDeflectingAod:(aodNumber-1))).*aodDirectionVectors{mthDeflectingAod}.*optimalBaseFreq; 
                xyDeflection = xyDeflection + offsetDueToMthAod;
            end
            xyDeflection = xyDeflection * lambda/V;
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
        
        f3diff = - V/lambda * xyDeflectionMm(1,:) * 1e-3 ./ (pairDeflectionRatio * sum(aodSpacing) + sum(aodSpacing(3:4)));
        f4diff = - V/lambda * xyDeflectionMm(2,:) * 1e-3 ./ (pairDeflectionRatio * sum(aodSpacing(2:4)) + aodSpacing(4));
        baseFreq = [ sum(aodSpacing(3:4))/sum(aodSpacing)*optimalBaseFreq - pairDeflectionRatio*f3diff,...
            aodSpacing(4)/sum(aodSpacing(2:4))*optimalBaseFreq - pairDeflectionRatio*f4diff, optimalBaseFreq+f3diff, optimalBaseFreq+f4diff];
    end
end