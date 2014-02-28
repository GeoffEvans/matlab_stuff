function [ prodEffOpt ] = aol_performance( microSecs, xMils, yMils, theta, phi, xDeflectMils, yDeflectMils, pairDeflectionRatio, optConstFreq, findFocusAndPlotRays, numToOptimise )

if size(theta,1) ~= size(phi,1) || size(theta,2) ~= size(phi,2) || length(xMils) ~= length(yMils) ||  length(xDeflectMils) ~= length(yDeflectMils)
    % size(theta) = size(phi) = [number of orientation perturbations, number of AODs]
    error('argument size mismatch');
end

[wavelengthVac, acPower, iPolAir, V] = AolConstants(); %TODO instantiate AODs with desired properties...
[numOfAods, aodAcDirectionVectors, zPlanesAod, linearChirps, constFreq, aodCentre] = AolDrive(theta, xDeflectMils, yDeflectMils, pairDeflectionRatio, optConstFreq); % TODO model the propagation through AOD thickness - adjust drive accordingly
[ prodEffOpt ] = AolPerformance( microSecs, xMils, yMils, theta, phi, findFocusAndPlotRays, numToOptimise );


% NOTE FOR MONDAY: this is broken because you tried to add in the different
% deflections. In doing so you removed the trace rays to determine the
% centre and now you are trying to find the centre by using a small angle
% approx. You need to get it plotting something that goes to a focus and
% then once you've done that check that a trace ray roughly matches up with
% the calculated centres. The number of arguments is getting to big so
% condense [xMils; yMils] and [xDefMil;yDefMil] to get from 11 to 9.

    function [wavelengthVac, acPower, iPolAir, V] = AolConstants()
        acPower = 1; % Watts
        iPolAir = [1; 1i]/sqrt(2);
        V = 613;
        wavelengthVac = 800e-9;
    end

    function [numOfAods, aodDirectionVectors, zPlanesAod, chirp, constFreq, aodCentre] = AolDrive(theta, xDeflection, yDeflection, pairDeflectionRatio, optConstFreq)
        numOfAods = size(theta,2);
        if numOfAods == 4
            [aodDirectionVectors, aodL, chirpFactor, constFreq] = Aod4(xDeflection, yDeflection, pairDeflectionRatio, optConstFreq);
        end
        if numOfAods == 2
            [aodDirectionVectors, aodL, chirpFactor] = Aod2();
        end
        if numOfAods == 1
            [aodDirectionVectors, aodL, chirpFactor] = Aod1();
        end
        chirp = V*V/wavelengthVac * chirpFactor;
        zPlanesAod = 1 + [0, cumsum(aodL)];
        
        aodCentre = zeros(2,numOfAods);
        for n = 2:numOfAods % calc centre of AODs using small angle approx.
            accumulator = 0;
            for m = 1:n
                offsetDueToMthAod = sum(aodL(m:(n-1))).*aodDirectionVectors{m}.*constFreq(m);
                accumulator = accumulator + offsetDueToMthAod;
            end
            aodCentre(:,n) = wavelengthVac/V * accumulator;
        end
        
        function [aodDirectionVectors, aodL, chirpFactor, constFreq] = Aod4(xDeflectMils, yDeflectMils, pairDeflectionRatio, optCentreFreq)
            %aodDirectionVectors = {[0;1], [1;0], -[0;1], -[1;0]};
            aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
            aodL = [5e-2, 5e-2, 5e-2, 5e0];
            l1 = aodL(1);
            l2 = aodL(2);
            l3 = aodL(3);
            l4 = aodL(4);
            chirpFactor = [ 1/(l1 + l2 + 2*l3 + 2*l4)...
                1/(l2 + l3 + 2*l4)...
                1/(2*(l3 + l4))...
                1/(2*l4)];      
            f3diff = - V/wavelengthVac * xDeflectMils ./ (pairDeflectionRatio * sum(aodL) + sum(aodL(3:4)));
            f4diff = - V/wavelengthVac * yDeflectMils ./ (pairDeflectionRatio * sum(aodL(2:4)) + aodL(4));
            constFreq = [ sum(aodL(3:4))/sum(aodL)*optCentreFreq - pairDeflectionRatio*f3diff; aodL(4)/sum(aodL(2:4))*optCentreFreq - pairDeflectionRatio*f4diff; optCentreFreq+f3diff; optCentreFreq+f4diff];
        end
        function [aodDirectionVectors, aodL, chirpFactor] = Aod2()
            aodDirectionVectors = {[1;0], -[1;0]};
            aodL = [5e-2, 5];
            l1 = aodL(1);
            l2 = aodL(2);
            chirpFactor = [ 1/(l1 + 2*l2), 1/(2*l2) ];
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
        PropagateThroughAods(k, theta, phi);
        zFocusModel = PropagateAfterAods(k,findFocusAndPlotRays);
        PlotRays(findFocusAndPlotRays,zFocusModel);
        [ prodEffOpt ] = AnalysePerformance(eff, numToOptimise, x(end-1,:), y(end-1,:), theta, phi);
        
        function [ prodEffOpt ] = AnalysePerformance(eff, numToOptimise, xFocus, yFocus, theta, phi)
            prodEffSingleTime = geomean(eff(1:numToOptimise,:),1); % average over all AODs we are interested in
            prodEffPerturbationAv = harmmean(reshape(prodEffSingleTime,numOfRaysPerPerturbation,numOfPerturbations),1); % average for each perturbation, want low efficiency to have big effect
            prodEffOpt = max(prodEffPerturbationAv);
            return % below not currently used but may be useful later / for debugging
            optIndices = prodEffTimeAv == max(prodEffTimeAv);
            sigmaX = std(reshape(xFocus,numOfRaysPerPerturbation,numOfPerturbations),1);
            sigmaY = std(reshape(yFocus,numOfRaysPerPerturbation,numOfPerturbations),1);
            avX = mean(reshape(xFocus,numOfRaysPerPerturbation,numOfPerturbations));
            avY = mean(reshape(yFocus,numOfRaysPerPerturbation,numOfPerturbations));
            maxFracAngleError = max(fractionalAngleErrorMax);
            thetaOpt = theta(optIndices,:) * 180/pi;
            phiOpt = phi(optIndices,:) * 180/pi;
        end
        
        function PropagateThroughAods(k, theta, phi)
            for a=1:numOfAods
                PropagateToAod(a, k, theta, phi);
                [k, eff(a,:)] = DeflectAtAod(a, k, theta, phi);
            end
            
            function PropagateToAod(aodNumber,k,theta,phi)
                % the vector from a point to anywhere on a plane dotted with the normal equals the distance to the plane from the point
                unitNormalToAod = GetUnitNormalToAod(aodNumber,theta,phi);
                vectorToCentreAod = StretchColumn([aodCentre(1,aodNumber); aodCentre(2,aodNumber); zPlanesAod(aodNumber)]) - [x(aodNumber,:); y(aodNumber,:); z(aodNumber,:)];
                x(aodNumber+1,:) = x(aodNumber,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(1,:);
                y(aodNumber+1,:) = y(aodNumber,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(2,:);
                z(aodNumber+1,:) = z(aodNumber,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(3,:);
                
                function normalUnitInLabFrame = GetUnitNormalToAod(aodNumber,theta,phi)
                    normalUnitInAodFrame = repmat([0; 0; 1],1,numOfRays);
                    normalUnitInLabFrame = TransformFrameAtNthAod(@PerturbedCrystalToLabFrame,normalUnitInAodFrame,aodNumber,theta,phi);
                end
            end
            
            function [kOut, eff] = DeflectAtAod(nthAod, kIn, theta, phi)
                Kn = aodAcDirectionVectors{nthAod};
                xFromCentre = x(nthAod,:) - aodCentre(1,nthAod);
                yFromCentre = y(nthAod,:) - aodCentre(2,nthAod);
                phase = StretchTimeArray(t) - ( xFromCentre*Kn(1) + yFromCentre*Kn(2) ) / V;
                localFreq = StretchBaseFrequencies(constFreq(nthAod,:)) + linearChirps(nthAod) * phase;
                
                kInR = TransformToPerturbedCrystalFrame(kIn,nthAod,theta,phi);
                [ kOutR, eff, ~ ] = aod3d.aod_propagator_vector( kInR, ones(1,numOfRays), StretchColumn(iPolAir), localFreq, StretchColumn(acPower) );
                kOut = TransformOutOfAdjustedCrystalFrame(kOutR,nthAod,theta,phi);
                
                modelAngle = acos( dot(kOut,kIn) ./ (mag(kOut).*mag(kIn)) );
                isotropicAngle = wavelengthVac * localFreq / V;
                fractionalAngleError = abs(isotropicAngle./modelAngle - 1);
                fractionalAngleErrorMax = fractionalAngleErrorMax + (fractionalAngleError > fractionalAngleErrorMax) .* (fractionalAngleError - fractionalAngleErrorMax);     
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
            zFocusModel = FindModelFocus(k,zFocusExpected,findFocus);
            [x(numOfAods+2,:),y(numOfAods+2,:),z(numOfAods+2,:)] = PropagateToNormalPlaneAfterLastAod(k, zFocusExpected); % expected focus
            [x(numOfAods+3,:),y(numOfAods+3,:),z(numOfAods+3,:)] = PropagateToNormalPlaneAfterLastAod(k, zFocusModel); % model focus
            [x(numOfAods+4,:),y(numOfAods+4,:),z(numOfAods+4,:)] = PropagateToNormalPlaneAfterLastAod(k, zFocusModel + 2); % past focus
            
            function focus = FindModelFocus(k,zFocusExpected,findFocus)
                if findFocus == true
                    kLocal = k;
                    focus = fminsearch(@MinFunc,zPlanesAod(end));
                else
                    focus = zFocusExpected;
                end
                function val = MinFunc(zVal)
                    [xTemp, yTemp, ~] = PropagateToNormalPlaneAfterLastAod(kLocal,zVal);
                    sigmaX = std(xTemp,1);
                    sigmaY = std(yTemp,1);
                    val = sigmaX .* sigmaY;
                end
            end
            
            function [xTemp,yTemp,zTemp] = PropagateToNormalPlaneAfterLastAod(k, zPosition)
                n = numOfAods+1;
                vectorToCentreOfPlane = -[x(n,:); y(n,:); z(n,:) - zPosition];
                xTemp = x(n,:) + vectorToCentreOfPlane(3,:) ./ k(3,:) .* k(1,:);
                yTemp = y(n,:) + vectorToCentreOfPlane(3,:) ./ k(3,:) .* k(2,:);
                zTemp = z(n,:) + vectorToCentreOfPlane(3,:);
            end
        end
        
        function PlotRays(plotRays,zFocusModel)
            if plotRays
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
            numOfDeflections = length(baseFreq);
            numOfPerturbations = size(theta,1);
            numOfRays = numOfPerturbations * numOfDeflections * numOfPositions * numOfTimes;
            numOfRaysPerTime = numOfRays / numOfTimes;
            numOfRaysPerPerturbation = numOfRays / numOfPerturbations;
        end
        function [t,x,y,z,eff,k] = InitialiseRayVars(microSecs,xMils,yMils)
            t = microSecs * 0;%1e-6;
            % the column structure of variables below should be [ Pertubrations{Deflections{Positions{Times}}} ]
            % the rows are input plane, AODs, expected focal plane, model focal plane, end plane - store these all for plotting
            x = zeros(numOfAods+4,numOfRays);
            x(1,:) = StretchPositionArray(xMils*1e-3);
            y = zeros(numOfAods+4,numOfRays);
            y(1,:) = StretchPositionArray(yMils*1e-3);
            z = zeros(numOfAods+4,numOfRays);
            eff = zeros(numOfAods,numOfRays);
            k = repmat([0;0;1]*2*pi/wavelengthVac,1,numOfRays); % input laser is orthogonal to AOD centre line
        end
        
        function [ reshaped ] = StretchTimeArray(array)
            reshaped = repmat(array,1,numOfRaysPerTime);
        end
        function [ reshaped ] = StretchPositionArray(array)
            reshaped = reshape(repmat(array,numOfTimes,numOfDeflections*numOfPerturbations),1,numOfRays);
        end
        function [ reshaped ] = StretchBaseFrequencies(array)
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