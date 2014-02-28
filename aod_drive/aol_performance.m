function [ prodEffOpt ] = aol_performance( microSecs, xMils, yMils, theta, phi, xDeflectMils, yDeflectMils, pairDeflectionRatio, findFocusAndPlotRays, numToOptimise )

if size(theta,1) ~= size(phi,1) || size(theta,2) ~= size(phi,2) || length(xMils) ~= length(yMils) ||  length(xDeflectMils) ~= length(yDeflectMils)
    % size(theta) = size(phi) = [number of orientation perturbations, number of AODs]
    error('argument size mismatch');
end
numOfAods = size(theta,2);

[wavelengthVac, acPower, iPolAir, V] = AolConstants(); %TODO instantiate AODs with desired properties...
[aodAcDirectionVectors, zPlanesAod, linearChirps, baseFreq] = AolDrive(numOfAods, xDeflectMils, yDeflectMils, pairDeflectionRatio);

[ ~,xOut,yOut ] = aol_performance_centred( 0, 0, 0, zeros(1,numOfAods), zeros(1,numOfAods), zeros(1,numOfAods), zeros(1,numOfAods), false, 1 ); % tracer ray to find AOD centres
[ prodEffOpt,~,~ ] = aol_performance_centred( microSecs, xMils, yMils, theta, phi, xOut(2:1+numOfAods), yOut(2:1+numOfAods), findFocusAndPlotRays, numToOptimise ); % do all rays

    function [wavelengthVac, acPower, iPolAir, V] = AolConstants()
        acPower = 1; % Watts
        iPolAir = [1; 1i]/sqrt(2);
        V = 613;
        wavelengthVac = 800e-9;
    end

    function [aodDirectionVectors, zPlanesAod, chirp, baseFreq] = AolDrive(number, xDeflection, yDeflection, pairDeflectionRatio)
        baseFreq = 30e6;
        if number == 4
            [aodDirectionVectors, zPlanesAod, chirp, baseFreq] = Aod4(xDeflection, yDeflection, pairDeflectionRatio);
        end
        if number == 2
            [aodDirectionVectors, zPlanesAod, chirp] = Aod2();
        end
        if number == 1
            [aodDirectionVectors, zPlanesAod, chirp] = Aod1();
        end
        
        function [aodDirectionVectors, zPlanesAod, chirp, baseFreq] = Aod4(xDeflection, yDeflection, pairDeflectionRatio)
            %aodDirectionVectors = {[0;1], [1;0], -[0;1], -[1;0]};
            aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
            aodL = [5e-2, 5e-2, 5e-2, 5e2];
            l1 = aodL(1);
            l2 = aodL(2);
            l3 = aodL(3);
            l4 = aodL(4);
            chirpFactor = [ 1/(l1 + l2 + 2*l3 + 2*l4)...
                1/(l2 + l3 + 2*l4)...
                1/(2*(l3 + l4))...
                1/(2*l4)];
            chirp = V*V/wavelengthVac * chirpFactor;
            zPlanesAod = 1 + [0, cumsum(aodL)];
            baseFreq = repmat(30e6,4,length(xDeflection)*length(pairDeflectionRatio));
        end
        function [aodDirectionVectors, zPlanesAod, chirp] = Aod2()
            aodDirectionVectors = {[1;0], -[1;0]};
            aodL = [5e-2, 5];
            l1 = aodL(1);
            l2 = aodL(2);
            chirpFactor = [ 1/(l1 + 2*l2), 1/(2*l2) ];
            chirp = V*V/wavelengthVac * chirpFactor;
            zPlanesAod = 1 + [0, cumsum(aodL)];
        end
        function [aodDirectionVectors, zPlanesAod, chirp] = Aod1()
            aodDirectionVectors = {[1;0]};
            aodL = 5;
            l1 = aodL(1);
            chirpFactor = 1/l1;
            chirp = V*V/wavelengthVac * chirpFactor;
            zPlanesAod = 1 + [0, cumsum(aodL)];
        end
    end

    function [ prodEffOpt,x,y ] = aol_performance_centred( microSecs, xMils, yMils, theta, phi, aodCentreX, aodCentreY, findFocusAndPlotRays, numToOptimise )
        fractionalAngleErrorMax = 0;
        [numOfTimes, numOfPositions, numOfDeflections, numOfPerturbations, numOfRays, numOfRaysPerTime, numOfRaysPerPerturbation] = UsefulNumbers(microSecs, xMils, theta, baseFreq);
        [t,x,y,z,eff,k] = InitialiseRayVars(microSecs,xMils,yMils);
        for a=1:numOfAods
            PropagateToAod(a, k, theta, phi);
            [k, eff(a,:)] = DeflectAtAod(a, k, theta, phi);
        end
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
        
        function [kOut, eff] = DeflectAtAod(nthAod, kIn, theta, phi)
            Kn = aodAcDirectionVectors{nthAod};
            phase = StretchTimeArray(t) - ( x(nthAod,:)*Kn(1) + y(nthAod,:)*Kn(2) ) / V;
            localFreq = StretchBaseFrequencies(baseFreq(nthAod,:)) + linearChirps(nthAod) * phase;
            
            kInR = TransformToPerturbedCrystalFrame(kIn,nthAod,theta,phi);
            [ kOutR, eff, ~ ] = aod3d.aod_propagator_vector( kInR, ones(1,numOfRays), StretchColumn(iPolAir), localFreq, StretchColumn(acPower) );
            kOut = TransformOutOfAdjustedCrystalFrame(kOutR,nthAod,theta,phi);
            
            modelAngle = acos( dot(kOut,kIn) ./ (mag(kOut).*mag(kIn)) );
            isotropicAngle = wavelengthVac * localFreq / V;
            fractionalAngleError = abs(isotropicAngle./modelAngle - 1);
            fractionalAngleErrorMax = fractionalAngleErrorMax + (fractionalAngleError > fractionalAngleErrorMax) .* (fractionalAngleError - fractionalAngleErrorMax);
            
            function kInR = TransformToPerturbedCrystalFrame(kIn,n,theta,phi)
                kInR = TransformFrameForNthAod(@LabToPerturbedCrystalFrame,kIn,n,theta,phi);
            end
            function kOut = TransformOutOfAdjustedCrystalFrame(kOutR,n,theta,phi)
                kOut = TransformFrameForNthAod(@PerturbedCrystalToLabFrame,kOutR,n,theta,phi);
            end
            function m = mag(v)
                m = dot(v,v);
                m = sqrt(m);
            end
        end
        
        function PropagateToAod(aodNumber,k,theta,phi)
            % the vector from a point to anywhere on a plane dotted with the normal equals the distance to the plane from the point
            unitNormalToAod = GetUnitNormalToAod(aodNumber,theta,phi);
            vectorToCentreAod = StretchColumn([aodCentreX(aodNumber); aodCentreY(aodNumber); zPlanesAod(aodNumber)]) - [x(aodNumber,:); y(aodNumber,:); z(aodNumber,:)];
            x(aodNumber+1,:) = x(aodNumber,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(1,:);
            y(aodNumber+1,:) = y(aodNumber,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(2,:);
            z(aodNumber+1,:) = z(aodNumber,:) + dot(vectorToCentreAod,unitNormalToAod) ./ dot(k,unitNormalToAod) .* k(3,:);
            
            function normalUnitInLabFrame = GetUnitNormalToAod(aodNumber,theta,phi)
                normalUnitInAodFrame = repmat([0; 0; 1],1,numOfRays);
                normalUnitInLabFrame = TransformFrameForNthAod(@PerturbedCrystalToLabFrame,normalUnitInAodFrame,aodNumber,theta,phi);
            end
        end
        
        function [xTemp,yTemp,zTemp] = PropagateToNormalPlaneAfterLastAod(k, zPosition)
            n = numOfAods+1;
            vectorToCentreOfPlane = -[x(n,:); y(n,:); z(n,:) - zPosition];
            xTemp = x(n,:) + vectorToCentreOfPlane(3,:) ./ k(3,:) .* k(1,:);
            yTemp = y(n,:) + vectorToCentreOfPlane(3,:) ./ k(3,:) .* k(2,:);
            zTemp = z(n,:) + vectorToCentreOfPlane(3,:);
        end
        
        function vOut = TransformFrameForNthAod(transformFun, vIn, n, theta, phi)
            vOut = zeros(3,numOfRays);
            for m = 1:numOfPerturbations
                phiAod = phi(m,n);
                thetaAod = theta(m,n);
                raysForMthPerm = (1:numOfRaysPerPerturbation)+(m-1)*numOfRaysPerPerturbation;
                vOut(:,raysForMthPerm) = transformFun(n,thetaAod,phiAod) * vIn(:,raysForMthPerm);
            end
        end
        
        function rotation = PerturbedCrystalToLabFrame(n,thetaAod,phiAod)
            rotation = transpose(LabToPerturbedCrystalFrame(n,thetaAod,phiAod));
        end
        function rotation = LabToPerturbedCrystalFrame(n,thetaAod,phiAod)
            rotation = UnperturbedToPerturbedCrystalFrame(thetaAod,phiAod) * LabToUnperturbedCrystalFrame(n);
        end
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

