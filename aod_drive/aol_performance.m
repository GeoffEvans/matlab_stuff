function [ prodEffOpt, x, y ] = aol_performance( microSecs, xMils, yMils, theta, phi, xDeflectMils, yDeflectMils, pairDeflectionRatio, optConstFreq, findFocusAndPlotRays, numToOptimise )

if size(theta,1) ~= size(phi,1) || size(theta,2) ~= size(phi,2) || length(xMils) ~= length(yMils) ||  size(xDeflectMils,2) ~= size(yDeflectMils,2) ||  size(xDeflectMils,1) ~= 1
    % size(theta) = size(phi) = [number of orientation perturbations, number of AODs]
    error('argument size mismatch');
end

[wavelengthVac, acPower, iPolAir, V] = AolConstants(); %TODO instantiate AODs with desired properties...
[numOfAods, aodAcDirectionVectors, zPlanesAod, linearChirps, constFreq, aodCentre] = AolDrive(theta, xDeflectMils, yDeflectMils, pairDeflectionRatio, optConstFreq); % TODO model the propagation through AOD thickness - adjust drive accordingly
[ prodEffOpt ] = AolPerformance( microSecs, xMils, yMils, theta, phi, findFocusAndPlotRays, numToOptimise );

    function [wavelengthVac, acPower, iPolAir, V] = AolConstants()
        acPower = 2; % Watts
        iPolAir = [1; 1i]/sqrt(2);
        V = 613;
        wavelengthVac = aod3d.opWavelenVac;
    end

    function [numOfAods, aodDirectionVectors, zPlanesAod, chirp, constFreq, aodCentre] = AolDrive(theta, xDeflection, yDeflection, pairDeflectionRatio, optConstFreq)
        correctionDistance = aod3d.L * ( 1 - 1/2.26 );
        numOfAods = size(theta,2);
        if numOfAods == 4
            [aodDirectionVectors, aodL, chirpFactor, constFreq] = Aod4(xDeflection, yDeflection, pairDeflectionRatio, optConstFreq);
        end
        if numOfAods == 2
            [aodDirectionVectors, aodL, chirpFactor, constFreq] = Aod2(xDeflection, pairDeflectionRatio, optConstFreq);
        end
        if numOfAods == 1
            [aodDirectionVectors, aodL, chirpFactor] = Aod1();
        end
        chirp = V*V/wavelengthVac * chirpFactor;
        zPlanesAod = 1 + [0, cumsum(aodL)];
        
        aodCentre = zeros(2,numOfAods+1);% calc centre of AODs and the focal point using small angle approx. 
        for n = 2:numOfAods+1 % start from 2 because offset on first AOD is zero
            accumulator = 0;
            for mAod = 1:n-1 
                offsetDueToMthAod = sum(aodL(mAod:(n-1))).*aodDirectionVectors{mAod}.*constFreq(mAod,1); % TODO adjust this to take const freq corresponding to zero deflection
                accumulator = accumulator + offsetDueToMthAod;
            end
            aodCentre(:,n) = wavelengthVac/V * accumulator;
        end
        
        function [aodDirectionVectors, aodL, chirpFactor, constFreq] = Aod4(xDeflectMils, yDeflectMils, pairDeflectionRatio, optConstFreq)
            %aodDirectionVectors = {[0;1], [1;0], -[0;1], -[1;0]};
            aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
            aodL = [5e-2, 5e-2, 5e-2, 1];
            
            l1 = aodL(1) - correctionDistance;
            l2 = aodL(2) - correctionDistance;
            l3 = aodL(3) - correctionDistance;
            l4 = aodL(4) - correctionDistance;
            chirpFactor = [ 1/(l1 + l2 + 2*l3 + 2*l4)...
                1/(l2 + l3 + 2*l4)...
                1/(2*(l3 + l4))...
                1/(2*l4)];      
            f3diff = - V/wavelengthVac * xDeflectMils * 1e-3 ./ (pairDeflectionRatio * sum(aodL) + sum(aodL(3:4)));
            f4diff = - V/wavelengthVac * yDeflectMils * 1e-3 ./ (pairDeflectionRatio * sum(aodL(2:4)) + aodL(4));
            constFreq = [ sum(aodL(3:4))/sum(aodL)*optConstFreq - pairDeflectionRatio*f3diff;...
                            aodL(4)/sum(aodL(2:4))*optConstFreq - pairDeflectionRatio*f4diff;...
                            optConstFreq+f3diff;...
                            optConstFreq+f4diff];
        end
        function [aodDirectionVectors, aodL, chirpFactor, constFreq] = Aod2(xDeflectMils, pairDeflectionRatio, optConstFreq)
            aodDirectionVectors = {[1;0], -[1;0]};
            aodL = [5e-2, 5];
            l1 = aodL(1) - correctionDistance;
            l2 = aodL(2) - correctionDistance;
            chirpFactor = [ 1/(l1 + 2*l2), 1/(2*l2) ];
            fDiff = - V/wavelengthVac * xDeflectMils * 1e-3 ./ (pairDeflectionRatio * sum(aodL) + aodL(2));
            constFreq = [ aodL(2)/sum(aodL)*optConstFreq - pairDeflectionRatio*fDiff; optConstFreq+fDiff];
        end
        function [aodDirectionVectors, aodL, chirpFactor] = Aod1()
            aodDirectionVectors = {[1;0]};
            aodL = 5;
            l1 = aodL(1);
            chirpFactor = 1/l1;
        end
    end

    function [ prodEffOpt ] = AolPerformance( microSecs, xMils, yMils, theta, phi, findFocusAndPlotRays, numToOptimise )
        fractionalAngleErrorMax = 0;
        [numOfTimes, numOfPositions, numOfDeflections, numOfPerturbations, numOfRays, numOfRaysPerTime, numOfRaysPerPerturbation] = UsefulNumbers(microSecs, xMils, theta, constFreq);
        [t,x,y,z,eff,k] = InitialiseRayVars(microSecs,xMils,yMils);
        k = PropagateThroughAods(k, theta, phi);
        zFocusModel = PropagateAfterAods(k,findFocusAndPlotRays);
        PlotRays(findFocusAndPlotRays,zFocusModel);
        [ prodEffOpt ] = AnalysePerformance(eff, numToOptimise, x(end-1,:), y(end-1,:), theta, phi);
        
        function [ prodEffDeflection ] = AnalysePerformance(eff, numToOptimise, xFocus, yFocus, theta, phi)
            prodEffSingleRay = geomean(eff(1:numToOptimise,:),1); % average over all AODs we are interested in
            prodEffDeflectionMat = reshape(prodEffSingleRay,numOfTimes*numOfPositions,numOfDeflections*numOfPerturbations); 
            prodEffDeflection = mean(prodEffDeflectionMat,1); % average for each deflection
            maxFracAngleError = max(fractionalAngleErrorMax);
        end
        
        function k = PropagateThroughAods(k, theta, phi)
            for aodNumber=1:numOfAods
                PropagateToAod(aodNumber, k, theta, phi);
                k = DeflectAtAod(aodNumber, k, theta, phi);
            end
            
            function PropagateToAod(aodNumber,k,theta,phi)                
                % the vector from a point to anywhere on a plane dotted with the normal equals the distance to the plane from the point
                unitNormalToAod = GetUnitNormalToAod(aodNumber,theta,phi);
                vectorToCentreAod = StretchColumn([aodCentre(1,aodNumber); aodCentre(2,aodNumber); zPlanesAod(aodNumber)]) - [x(aodNumber*2-1,:); y(aodNumber*2-1,:); z(aodNumber*2-1,:)];
                aodEntrance = IndexForAodFront(aodNumber);
                x(aodEntrance,:) = x(aodEntrance-1,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(1,:);
                y(aodEntrance,:) = y(aodEntrance-1,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(2,:);
                z(aodEntrance,:) = z(aodEntrance-1,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(3,:);
                
                function normalUnitInLabFrame = GetUnitNormalToAod(aodNumber,theta,phi)
                    normalUnitInAodFrame = repmat([0; 0; 1],1,numOfRays);
                    normalUnitInLabFrame = TransformFrameAtNthAod(@PerturbedCrystalToLabFrame,normalUnitInAodFrame,aodNumber,theta,phi);
                end
            end
            
            function kOut = DeflectAtAod(nthAod, kIn, theta, phi)
                localFreq = FindLocalPhase(nthAod);
                [thetaBragg, phiBragg] = CalculateBraggRotationAngles(kIn, aodAcDirectionVectors{nthAod}, localFreq); % for comparison with optimisations
                kInR = TransformToPerturbedCrystalFrame(kIn,nthAod,theta,phi);
                [ dispInCrystal, kOutR, eff(nthAod,:), ~ ] = aod3d.aod_propagator_vector( kInR, ones(1,numOfRays), StretchColumn(iPolAir), localFreq, StretchColumn(acPower) );
                kOut = TransformOutOfAdjustedCrystalFrame(kOutR,nthAod,theta,phi);
                dispInLab = TransformOutOfAdjustedCrystalFrame(dispInCrystal,nthAod,theta,phi);
                CalculateRayPositionAtAodExit(dispInLab, nthAod);
                CompareModelToSimpleAngles(kIn, kOut, localFreq);
                
                function localFreq = FindLocalPhase(nthAod)
                    Kn = aodAcDirectionVectors{nthAod};
                    aodEntrance = IndexForAodFront(nthAod);
                    xFromCentre = x(aodEntrance,:) - aodCentre(1,nthAod);
                    yFromCentre = y(aodEntrance,:) - aodCentre(2,nthAod);
                    phase = StretchTimeArray(t) - ( xFromCentre*Kn(1) + yFromCentre*Kn(2) ) / V;
                    localFreq = StretchDeflections(constFreq(nthAod,:)) + linearChirps(nthAod) * phase;
                end
                
                function CalculateRayPositionAtAodExit(dispInLab, nthAod)
                    aodEntrance = IndexForAodFront(nthAod);
                    x(aodEntrance+1,:) = x(aodEntrance,:) + dispInLab(1,:);
                    y(aodEntrance+1,:) = y(aodEntrance,:) + dispInLab(2,:);
                    z(aodEntrance+1,:) = z(aodEntrance,:) + dispInLab(3,:);
                end
                
                function CompareModelToSimpleAngles(kIn, kOut, localFreq)
                    modelAngle = acos( dot(kOut,kIn) ./ (mag(kOut).*mag(kIn)) );
                    isotropicAngle = wavelengthVac * localFreq / V;
                    fractionalAngleError = abs((isotropicAngle - modelAngle)./modelAngle);
                    fractionalAngleErrorMax = fractionalAngleErrorMax + (fractionalAngleError > fractionalAngleErrorMax) .* (fractionalAngleError - fractionalAngleErrorMax);  
                end
                
                function [theta,phi] = CalculateBraggRotationAngles(k,unitK2d,freq)
                    V = 613;
                    unitK = [unitK2d;0];
                    unitKarr = repmat(unitK,1,numOfRays);
                    zxK = cross([0;0;1],unitK);
                    zxKarr = repmat(zxK,1,numOfRays);
                    k2d = k - repmat(dot(k,zxKarr),3,1).*zxKarr;
                    k2dMag = mag(k2d);
                    theta = abs( acos(k(3)./k2dMag).*2.*( (dot(unitKarr,k)<0) - 0.5 ) - asin(aod3d.opWavelenVac*freq/2/V) );
                    phi = acos( dot(k,unitKarr) ./ mag(k) ) .* 2.*( (dot(zxKarr,k)>0) - 0.5 );
                end
            end
            
            function kInR = TransformToPerturbedCrystalFrame(kIn,n,theta,phi)
                kInR = TransformFrameAtNthAod(@LabToPerturbedCrystalFrame,kIn,n,theta,phi);
            end
            function kOut = TransformOutOfAdjustedCrystalFrame(kOutR,n,theta,phi)
                kOut = TransformFrameAtNthAod(@PerturbedCrystalToLabFrame,kOutR,n,theta,phi);
            end
            function rotation = PerturbedCrystalToLabFrame(n,thetaAod,phiAod)
                rotation = transpose(LabToPerturbedCrystalFrame(n,thetaAod,phiAod));
            end
            function rotation = LabToPerturbedCrystalFrame(n,thetaAod,phiAod)
                rotation = UnperturbedToPerturbedCrystalFrame(thetaAod,phiAod) * LabToUnperturbedCrystalFrame(n);
                
                function rotation = LabToUnperturbedCrystalFrame(n)
                    % change reference frame to have acoustic vector down <110> and z remains <001>
                    Kn = aodAcDirectionVectors{n};
                    Kx = Kn(1);
                    Ky = Kn(2);
                    rotation = [1 -1 0; 1 1 0; 0 0 sqrt(2)]/sqrt(2)*[Kx, Ky 0; -Ky, Kx 0; 0 0 1];
                end
                
                function rotation = UnperturbedToPerturbedCrystalFrame(thetaAod,phiAod)
                    % rotate the crystal by theta about the axis at angle phi to the acoustic vector in the xy plane (no rotation about z), change to this new frame
                    rotation = RotationPhi(phiAod-pi/4) * RotationTheta(thetaAod)' * RotationPhi(phiAod-pi/4)';
                    
                    function mat = RotationPhi(phi)
                        % from x towards y about z
                        mat = [cos(phi) -sin(phi) 0;...
                            sin(phi) cos(phi) 0;...
                            0           0     1];
                    end
                    function mat = RotationTheta(theta)
                        % from x to z about y
                        mat = [cos(theta) 0 -sin(theta);...
                            0       1       0;...
                            sin(theta)  0   cos(theta)];
                    end
                end
            end
        end
        
        function [zFocusModel] = PropagateAfterAods(k, findFocus)
            zFocusExpected = zPlanesAod(end);
            lastAodExitIndex = IndexForAodFront(numOfAods) + 1;
            zFocusModel = FindModelFocus(k,zFocusExpected,findFocus, lastAodExitIndex);
            [x(lastAodExitIndex+1,:),y(lastAodExitIndex+1,:),z(lastAodExitIndex+1,:)] = PropagateToNormalPlaneAfterLastAod(k, zFocusExpected, lastAodExitIndex); % expected focus
            [x(lastAodExitIndex+2,:),y(lastAodExitIndex+2,:),z(lastAodExitIndex+2,:)] = PropagateToNormalPlaneAfterLastAod(k, zFocusModel, lastAodExitIndex); % model focus
            [x(lastAodExitIndex+3,:),y(lastAodExitIndex+3,:),z(lastAodExitIndex+3,:)] = PropagateToNormalPlaneAfterLastAod(k, zFocusModel + zFocusExpected/4, lastAodExitIndex); % past focus
            
            function focus = FindModelFocus(k,zFocusExpected,findFocus, lastAodExitIndex)
                if findFocus == true
                    kLocal = k;
                    focus = fminsearch(@MinFunc,zPlanesAod(end));
                else
                    focus = zFocusExpected;
                end
                function val = MinFunc(zVal)
                    [xTemp, yTemp, ~] = PropagateToNormalPlaneAfterLastAod(kLocal,zVal,lastAodExitIndex);
                    xDeflectionCol = reshape(permute(reshape(xTemp,numOfTimes*numOfPositions,numOfDeflections,numOfPerturbations),[1 3 2]),numOfTimes*numOfPositions*numOfPerturbations,numOfDeflections);
                    yDeflectionCol = reshape(permute(reshape(yTemp,numOfTimes*numOfPositions,numOfDeflections,numOfPerturbations),[1 3 2]),numOfTimes*numOfPositions*numOfPerturbations,numOfDeflections);
                    sigmaX = std(xDeflectionCol,1);
                    sigmaY = std(yDeflectionCol,1);
                    val = prod(sigmaX) .* prod(sigmaY);
                end
            end
            
            function [xPos,yPos,zPos] = PropagateToNormalPlaneAfterLastAod(k, zPosition, lastAodExitIndex)
                vectorToCentreOfPlane = -[x(lastAodExitIndex,:); y(lastAodExitIndex,:); z(lastAodExitIndex,:) - zPosition];
                xPos = x(lastAodExitIndex,:) + vectorToCentreOfPlane(3,:) ./ k(3,:) .* k(1,:);
                yPos = y(lastAodExitIndex,:) + vectorToCentreOfPlane(3,:) ./ k(3,:) .* k(2,:);
                zPos = z(lastAodExitIndex,:) + vectorToCentreOfPlane(3,:);
            end
        end
        
        function PlotRays(plotRays,zFocusModel)
            if plotRays
                figure()
                hold on;
                for m=1:numOfAods+1
                    zValAsArray = repmat(zPlanesAod(m),1,4);
                    fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
                end
                zValAsArray = repmat(zFocusModel,1,4);
                fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
                alpha(0.1)
                for q = 1:numOfPerturbations
                    indicesForQthPerturbation = (1:numOfRaysPerPerturbation)+(q-1)*numOfRaysPerPerturbation;
                    plot3(x(:,indicesForQthPerturbation),y(:,indicesForQthPerturbation),z(:,indicesForQthPerturbation),'r');
                end
                grid on;
                grid minor;
                axis square;
                xlabel('x')
                ylabel('y')
                zlabel('z')
                hold off;
            end
        end
        
        function vOut = TransformFrameAtNthAod(matrixFromAnglesFunc, vIn, n, theta, phi)
            vOut = zeros(3,numOfRays);
            for m = 1:numOfPerturbations
                phiAod = phi(m,n);
                thetaAod = theta(m,n);
                raysForMthPerm = (1:numOfRaysPerPerturbation)+(m-1)*numOfRaysPerPerturbation;
                vOut(:,raysForMthPerm) = matrixFromAnglesFunc(n,thetaAod,phiAod) * vIn(:,raysForMthPerm);
            end
        end
           
        function [numOfTimes, numOfPositions, numOfDeflections, numOfPerturbations, numOfRays, numOfRaysPerTime, numOfRaysPerPerturbation] = UsefulNumbers(microSecs, xMils, theta, baseFreq)
            numOfTimes = length(microSecs);
            numOfPositions = length(xMils);
            numOfDeflections = size(baseFreq,2);
            numOfPerturbations = size(theta,1);
            numOfRays = numOfPerturbations * numOfDeflections * numOfPositions * numOfTimes;
            numOfRaysPerTime = numOfRays / numOfTimes;
            numOfRaysPerPerturbation = numOfRays / numOfPerturbations;
        end
        function [t,x,y,z,eff,k] = InitialiseRayVars(microSecs,xMils,yMils)
            t = microSecs * 0;%1e-6;
            % the column structure of variables below should be [ Pertubrations{Deflections{Positions{Times}}} ]
            % the rows are input plane, [AOD front AOD back] * 4, expected focal plane, model focal plane, end plane - store these all for plotting
            x = zeros(numOfAods*2+4,numOfRays);
            x(1,:) = StretchPositionArray(xMils*1e-3);
            y = zeros(numOfAods*2+4,numOfRays);
            y(1,:) = StretchPositionArray(yMils*1e-3);
            z = zeros(numOfAods*2+4,numOfRays);
            eff = zeros(numOfAods,numOfRays);
            k = repmat([0;0;1]*2*pi/wavelengthVac,1,numOfRays); % input laser is orthogonal to AOD centre line
        end
        function idx = IndexForAodFront(n)
            idx = n * 2;
        end
        
        function [ reshaped ] = StretchTimeArray(array)
            reshaped = repmat(array,1,numOfRaysPerTime);
        end
        function [ reshaped ] = StretchPositionArray(array)
            reshaped = reshape(repmat(array,numOfTimes,numOfDeflections*numOfPerturbations),1,numOfRays);
        end
        function [ reshaped ] = StretchDeflections(array)
            reshaped = reshape(repmat(array,numOfTimes*numOfPositions,numOfPerturbations),1,numOfRays);
        end
        function [ reshaped ] = StretchColumn(vector)
            reshaped = repmat(vector,1,numOfRays);
        end
    end
end

function m = mag(v)
m = dot(v,v);
m = sqrt(m);
end