function [ rayBundle ] = aol_model( microSecs, xyInputMm, thetaPhiAodPerturbations, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed )

    numOfAods = size(thetaPhiAodPerturbations,2);
    
    [aodDirectionVectors, aodCentres, aolDrives] = aol_drive_parameters(numOfAods, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed);
    
    rayBundle = aol_ray_bundle(microSecs,xyInputMm,aolDrives,thetaPhiPerturbs);    
    PropagateThroughAods(rayBundle);    
    trace_rays_after_aol(rayBundle, zFocusExpected)
    
    function rayBundle = PropagateThroughAods(numOfAods, rayBundle)
        
        for aodNumber=1:numOfAods
            PropagateToAod(aodNumber, rayBundle);
            DeflectAtAod(aodNumber, rayBundle);
        end
        
        function PropagateToAod(aodNumber,rayBundle)
            unitNormalToAod = GetUnitNormalToAod(aodNumber,rayBundle);
            xyzAod = propagate_ray_to_plane(rayBundle.GetXyzForAodBack(aodNumber),rayBundle.k,unitNormalToAod,aodCentres(:,aodNumber));
            rayBundle.SetXyzForAodFront(aodNumber,xyzAod);
            
            function normalUnitInLabFrame = GetUnitNormalToAod(aodNumber,rayBundle)
                normalUnitInAodFrame = repmat([0; 0; 1],rayBundle.numOfPerturbations);
                normalUnitInLabFrame = TransformFrameAtNthAod(@PerturbedCrystalToLabFrame,normalUnitInAodFrame,aodNumber,rayBundle.perturbsTheta,rayBundle.perturbsPhi);
                % Need to scale up from number of perturbations to number of rays so reshape:
                normalUnitInLabFrame = repmat(normalUnitInLabFrame,rayBundle.numOfRaysPerPerturbations,1);
                normalUnitInLabFrame = reshape(normalUnitInLabFrame,3,rayBundle.numOfRaysPerPerturbations,rayBundle.numOfPerturbations);
            end
        end
        
        function DeflectAtAod(nthAod, rayBundle)
            theta = rayBundle.perturbsTheta;
            phi = rayBundle.perturbsPhi;
            
            localFreq = FindLocalPhase(nthAod, rayBundle, theta, phi);

            kInCf = TransformToPerturbedCrystalFrame(rayBundle.k, nthAod, theta, phi);
            [ displacementInCrystalCf, kOutCf, rayBundle.eff(nthAod,:) ] = AodModel( kInCf, localFreq );
            rayBundle.k = TransformOutOfPerturbedCrystalFrame(kOutCf, nthAod, theta, phi);
            
            displacementInCrystal = TransformOutOfPerturbedCrystalFrame(displacementInCrystalCf,nthAod,theta,phi);
            rayBundle.SetXyzNthAodBack(nthAod, xyzIn + displacementInCrystal);
            
            function localFreq = FindLocalPhase(nthAod, rayBundle, theta, phi)
                xyzIn = rayBundle.GetXzyForAodFront(nthAod);
                xyzInCf = TransformToPerturbedCrystalFrame(xyzIn,nthAod,theta,phi);
                xFromCentreCf = xyzInCf(1,:) - aodCentres(1,nthAod);
                yFromCentreCf = xyzInCf(2,:) - aodCentres(2,nthAod);
                
                acousticDirectionCf = [1 1 0];
                phase = time - ( xFromCentreCf*acousticDirectionCf(1) + yFromCentreCf*acousticDirectionCf(2) ) / V;
                localFreq = rayBundle.drives.baseFreq(nthAod) + rayBundle.drives.chirp(nthAod) * phase;
            end
        end
        
        function kInR = TransformToPerturbedCrystalFrame(kIn,n,theta,phi)
            kInR = TransformFrameAtNthAod(@LabToPerturbedCrystalFrame,kIn,n,theta,phi);
        end
        function kOut = TransformOutOfPerturbedCrystalFrame(kOutR,n,theta,phi)
            kOut = TransformFrameAtNthAod(@PerturbedCrystalToLabFrame,kOutR,n,theta,phi);
            function rotation = PerturbedCrystalToLabFrame(n,thetaAod,phiAod)
                rotation = transpose(LabToPerturbedCrystalFrame(n,thetaAod,phiAod));
            end
        end
        function rotation = LabToPerturbedCrystalFrame(n,thetaAod,phiAod)
            rotation = UnperturbedToPerturbedCrystalFrame(thetaAod,phiAod) * LabToUnperturbedCrystalFrame(n);
            
            function rotation = LabToUnperturbedCrystalFrame(n)
                % change reference frame to have acoustic vector down <110> and z remains <001>
                Kn = aodDirectionVectors{n};
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
        function vectorsOut = TransformFrameAtNthAod(MapAngleToMatrix, vectorsIn, nthAod, theta, phi)
            vectorsOut = zeros(3,numOfRaysPerPerturbation,numOfPerturbations);
            for m = 1:numOfPerturbations
                phiAod = phi(m,nthAod);
                thetaAod = theta(m,nthAod);
                raysForMthPerm = (1:numOfRaysPerPerturbation)+(m-1)*numOfRaysPerPerturbation;
                vectorsOut(:,raysForMthPerm) = MapAngleToMatrix(nthAod,thetaAod,phiAod) * vectorsIn(:,raysForMthPerm);
            end
        end
    end
end

function [ dispInCrystal, kOutR, eff ] = AodModel( kInR, localFreq )
    numOfRays = length(localFreq);
    acPower = 2; % Watts
    iPolAir = [1; 1i]/sqrt(2); % Circular
    [ dispInCrystal, kOutR, eff ] = aod3d.aod_propagator_vector( kInR, ones(1,numOfRays), repmat(iPolAir,1,numOfRays), localFreq, repmat(acPower,1,numOfRays) );
end

