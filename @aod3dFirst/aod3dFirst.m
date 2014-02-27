classdef aod3dFirst
    
    properties (Constant)
        L = 3.6e-3; % transducer length
        opWavelenVac = 800e-9;
    end
    
    methods (Static)        % Externally defined
        [ acTheta, dTheta, dPhi, acInverseWavelenOpt, nOrdOpt, nExt ] =...
            match_phase( iTheta, iPhi, acFreq )
        [eff,dTheta,dPhi] = rescattered_efficiency( iTheta, iPhi, acFreq )
        [ dThetaAir, dPhiAir, dIntensity, dPolarisation ] = ...
            aod_propagator( iThetaAir, iPhiAir, acFreq, iIntensity, iPolarisation )
        plot_efficiency_surface_qwp()
        plot_diffraction_angles_surface()
        plot_efficiency_surface_2aod()

        function [ extendedX, extendedY ] = SetCrossProduct( X, Y )
            [x, y] = meshgrid(X,Y);
            extendedX = x(:)';
            extendedY = y(:)';
        end
        function [ theta, phi] = OpticalRotationAngles(opticRotation, rotatedTheta)
            % optic rotation angle defined TOWARDS [-110] (y');
            % Use z, [110], [-110]
            v110 = sin(rotatedTheta);
            v001 = cos(rotatedTheta);
            vN10 = 0; % for completeness and method transparency
            % rotate in z -- [-110] plane
            v110oR = v110;
            v001oR = v001 .* cos(opticRotation);
            vN10oR = v001 .* sin(opticRotation);
            % [110] unchanged
            % rotate to x,y,z
            v100out = (v110oR - vN10oR)./sqrt(2);
            v010out = (v110oR + vN10oR)./sqrt(2);
            v001out = v001oR;
            % calculate theta, phi
            [~, theta, phi] = get_angles_from_vector([v100out; v010out; v001out]);
        end
        function [ theta, phi] = ConvertAnglesBetweenFrames(thetaIn, phiIn, thetaZ, phiZ, thetaX, phiX)
            % anglesX,Z indicate the orientations of the new axes in the old coords
            vX = get_vector_from_angles(1,thetaX,phiX);
            vZ = get_vector_from_angles(1,thetaZ,phiZ);
            vY = cross(vZ,vX);
            len = get_angles_from_vector(vY);
            if abs(len-1) > 1e-6
                error('Your X and Z axes arne`t orthogonal! Aborting...');
            end
            vIn = get_vector_from_angles(1,thetaIn,phiIn);
            % Now find the components of vIn in the basis of {vX,vY,vZ}
            vOut = vIn; % assign size
            vOut(1,:) = dot(vIn,vX);
            vOut(2,:) = dot(vIn,vY);
            vOut(3,:) = dot(vIn,vZ);
            [ ~, theta, phi ] = get_angles_from_vector(vOut);
        end
        function [ theta, phi ] = MakeAllThetaPositive(thetaIn, phiIn)
            negativeTheta = thetaIn < 0;
            multiplier = negativeTheta * -2 + 1;
            theta = multiplier .* thetaIn;
            phi = mod(phiIn + pi*negativeTheta, 2*pi);
            phiOverPi = phi > pi;
            phi = phi - (2*pi * phiOverPi);
        end
    end
    
end

