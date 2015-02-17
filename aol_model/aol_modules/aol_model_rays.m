function [ rb ] = aol_model_rays( microSecs, xyInputMm, aolPerturbations, driveParams, transducerWidths )
% Translates a bunch of parameters into rays and propagates them through an
% AOL, calculating the focal plane if parameters correspond to focusing
% at one point with a real AOL
   
    [aodDirectionVectors, aodCentres, zFocusPredicted, aolDrives] = calculate_aol_drive(aolPerturbations.numOfAods, driveParams, 1);
    
    rb = aol_ray_bundle(microSecs,xyInputMm,aolDrives,aodCentres,zFocusPredicted,aolPerturbations,driveParams.opWavelenVac);    
    PropagateThroughAods(rb);
    
    isPointingModeAndSingleBundle = (rb.numOfPerturbations * rb.numOfDrives == 1) && (abs(sum(driveParams.xyScanSpeed)) == 0);
    trace_rays_after_aol(rb, isPointingModeAndSingleBundle);
    
    function rb = PropagateThroughAods(rb)
        
        for aodNumber = 1:rb.numOfAods
            PropagateToAod(aodNumber, rb);
            DeflectAtAod(aodNumber, rb);
        end
        
        function PropagateToAod(aodNumber,rb)
            unitNormalToAod = GetUnitNormalToAod(aodNumber,rb);
            aodCentresMatrix = rb.aodCentres{aodNumber};
            xyzAod = propagate_ray_to_plane(rb.GetXyzNthAodBack(aodNumber-1),rb.k,unitNormalToAod,aodCentresMatrix);
            rb.SetXyzNthAodFront(aodNumber,xyzAod);
            
            function normalUnitInLabFrame = GetUnitNormalToAod(aodNumber,rb)
                normalUnitInAodFrameForEachPerturbation = repmat([0; 0; 1],[1,1,rb.numOfPerturbations]);
                normalUnitInLabFrameForEachPerturbation = TransformOutOfPerturbedCrystalFrame(normalUnitInAodFrameForEachPerturbation,aodNumber,rb);
                normalUnitInLabFrame = stretch(normalUnitInLabFrameForEachPerturbation,rb.numOfRaysPerPerturbation);
            end
        end
        
        function DeflectAtAod(nthAod, rb)
            localFreq = FindLocalFreq(nthAod, rb);

            kInCf = TransformToPerturbedCrystalFrame(rb.k, nthAod, rb);
            [ displacementInCrystalCf, kOutCf, rb.eff(nthAod,:,:) ] = AodModel( kInCf, localFreq, transducerWidths(nthAod), rb.opWavelenVac );
            rb.k = TransformOutOfPerturbedCrystalFrame(kOutCf, nthAod, rb);
            
            displacementInCrystal = TransformOutOfPerturbedCrystalFrame(displacementInCrystalCf,nthAod,rb);
            rb.SetXyzNthAodBack(nthAod, displacementInCrystal);
            
            function localFreq = FindLocalFreq(nthAod, rb)
                xyzIn = rb.GetXyzNthAodFront(nthAod);
                xyzFromCentre = xyzIn - rb.aodCentres{nthAod};
                xyzFromCentreCf = TransformToPerturbedCrystalFrame(xyzFromCentre,nthAod,rb);
                
                V = teo2.find_v_ac_min(pi/2,pi/4);
                acousticDirectionCf = [1; 1]/sqrt(2);
                
                phase = rb.t - ( xyzFromCentreCf(1,:,:)*acousticDirectionCf(1)...
                    + xyzFromCentreCf(2,:,:)*acousticDirectionCf(2) ) / V;
                
                localFreq = rb.BaseFreqForNthAod(nthAod) + rb.ChirpForNthAod(nthAod) .* phase;
            end
        end
        
        function kInR = TransformToPerturbedCrystalFrame(kIn,n,rb)
            kInR = rb.ApplyPerturbationMatricesToVectors(@LabToPerturbedCrystalFrame,kIn,n);
        end
        function kOut = TransformOutOfPerturbedCrystalFrame(kOutR,n,rb)
            kOut = rb.ApplyPerturbationMatricesToVectors(@PerturbedCrystalToLabFrame,kOutR,n);
            if ~isreal(sum(kOut))
                error(['imaginary angles have appeared propagating out of AOD ' num2str(n)])
            end
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

function [ dispInCrystal, kOut, eff ] = AodModel( kIn, localFreq, transducerWidth, opWavelenVac )
    numOfRays = numel(localFreq);
    acPower = 0.9; % Watts
    iPolAir = [1; 1i]/sqrt(2); % Circular
    localFreqSquashed = reshape(localFreq,1,numOfRays);
    kInSquashed = reshape(kIn,3,numOfRays);
    [ dispInCrystalSquashed, kOutSquashed, effSquashed ] = aod3d.aod_propagator_vector( kInSquashed, ones(1,numOfRays), repmat(iPolAir,1,numOfRays), localFreqSquashed, repmat(acPower,1,numOfRays), transducerWidth, opWavelenVac );
    dispInCrystal = reshape(dispInCrystalSquashed, size(kIn));
    kOut = reshape(kOutSquashed, size(kIn));
    eff = reshape(effSquashed, size(localFreq));
end

