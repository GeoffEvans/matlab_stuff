function [ rayBundle ] = aol_model( microSecs, xyInputMm, thetaPhiAodPerturbations, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed )
    
    numOfAods = size(thetaPhiAodPerturbations{1},2);
    
    [aodDirectionVectors, aodCentres, zFocusPredicted, aolDrives] = calculate_aol_drive(numOfAods, optimalBaseFreq, xyDeflectionMm, pairDeflectionRatio, scanSpeed);
    
    rayBundle = aol_ray_bundle(microSecs,xyInputMm,aolDrives,aodCentres,zFocusPredicted,thetaPhiAodPerturbations);    
    PropagateThroughAods(rayBundle);
    
    isPointingModeAndSingleBundle = (rayBundle.numOfPerturbations * rayBundle.numOfDrives == 1) && (scanSpeed == 0);
    rayBundle.zFocusModel = trace_rays_after_aol(rayBundle, zFocusExpected, isPointingModeAndSingleBundle);
    
    function rayBundle = PropagateThroughAods(rayBundle)
        
        for aodNumber = 1:rayBundle.numOfAods
            PropagateToAod(aodNumber, rayBundle);
            DeflectAtAod(aodNumber, rayBundle);
        end
        
        function PropagateToAod(aodNumber,rayBundle)
            unitNormalToAod = GetUnitNormalToAod(aodNumber,rayBundle);
            aodCentresMatrix = repmat(rayBundle.aodCentres(:,aodNumber),[1,rayBundle.numOfRaysPerPerturbation,rayBundle.numOfPerturbations]);
            xyzAod = propagate_ray_to_plane(rayBundle.GetXyzNthAodBack(aodNumber-1),rayBundle.k,unitNormalToAod,aodCentresMatrix);
            rayBundle.SetXyzNthAodFront(aodNumber,xyzAod);
            
            function normalUnitInLabFrame = GetUnitNormalToAod(aodNumber,rayBundle)
                normalUnitInAodFrame = repmat([0; 0; 1],rayBundle.numOfPerturbations);
                normalUnitInLabFrame = TransformOutOfPerturbedCrystalFrame(normalUnitInAodFrame,aodNumber,rayBundle);
                normalUnitInLabFrame = repmat(stretch(normalUnitInLabFrame,rayBundle.numOfRaysPerPerturbation),[1,1,rayBundle.numOfPerturbations]);
            end
        end
        
        function DeflectAtAod(nthAod, rayBundle)
            theta = rayBundle.perturbsTheta;
            phi = rayBundle.perturbsPhi;
            
            localFreq = FindLocalPhase(nthAod, rayBundle);

            kInCf = TransformToPerturbedCrystalFrame(rayBundle.k, nthAod, theta, phi);
            [ displacementInCrystalCf, kOutCf, rayBundle.eff(nthAod,:) ] = AodModel( kInCf, localFreq );
            rayBundle.k = TransformOutOfPerturbedCrystalFrame(kOutCf, nthAod, theta, phi);
            
            displacementInCrystal = TransformOutOfPerturbedCrystalFrame(displacementInCrystalCf,nthAod,theta,phi);
            rayBundle.SetXyzNthAodBack(nthAod, displacementInCrystal);
            
            function localFreq = FindLocalPhase(nthAod, rayBundle)
                xyzIn = rayBundle.GetXyzNthAodFront(nthAod);
                xyzInCf = TransformToPerturbedCrystalFrame(xyzIn,nthAod,rayBundle);
                xFromCentreCf = xyzInCf(1,:) - rayBundle.aodCentres(1,nthAod);
                yFromCentreCf = xyzInCf(2,:) - rayBundle.aodCentres(2,nthAod);
                
                V = teo2.find_v_ac_min(pi/2,pi/4);
                acousticDirectionCf = [1 1 0];
                
                phase = rayBundle.t - ( xFromCentreCf*acousticDirectionCf(1) + yFromCentreCf*acousticDirectionCf(2) ) / V;
                
                localFreq = rayBundle.BaseFreqForNthAod(nthAod) + rayBundle.ChirpForNthAod(nthAod) * phase;
            end
        end
        
        function kInR = TransformToPerturbedCrystalFrame(kIn,n,rayBundle)
            kInR = rayBundle.ApplyPerturbationMatricesToVectors(@LabToPerturbedCrystalFrame,kIn,n);
        end
        function kOut = TransformOutOfPerturbedCrystalFrame(kOutR,n,rayBundle)
            kOut = rayBundle.ApplyPerturbationMatricesToVectors(@PerturbedCrystalToLabFrame,kOutR,n);
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
    end
end

function [ dispInCrystal, kOutR, eff ] = AodModel( kInR, localFreq )
    numOfRays = length(localFreq);
    acPower = 2; % Watts
    iPolAir = [1; 1i]/sqrt(2); % Circular
    [ dispInCrystal, kOutR, eff ] = aod3d.aod_propagator_vector( kInR, ones(1,numOfRays), repmat(iPolAir,1,numOfRays), localFreq, repmat(acPower,1,numOfRays) );
end

