function [ acTheta, dTheta, dPhi, acInvWavelen, nOrd, nExt ] =...
    match_phase( iTheta, iPhi, acFreq )
% Finds the diffracted and acoustic beam angles for arrays of incident
% angles and acoustic frequencies.

% Error checking: all inputs same length arrays
% 0 < Theta < pi

[acTheta, dTheta, dPhi, acInvWavelen, nOrd] = SolveVectorTriangle(acFreq, iTheta, iPhi);

    function [acTheta, dTheta, dPhi, acInvWavelen, nOrd] = SolveVectorTriangle(acFreq, iTheta, iPhi)
        opWavelenVac = aod3dFirst.opWavelenVac;
        [ ~, nExt, ~, ~ ] = teo2.find_n_op( iTheta );
        k_i_norm = 2 * pi * nExt / opWavelenVac;
        k_i = get_vector_from_angles(k_i_norm, iTheta, iPhi);
        
        options = optimset('Algorithm','trust-region-reflective','display','off',...
            'JacobPattern',speye(length(iTheta)));
        acThetaStart = pi/2 + 0*iTheta;
        acTheta = fsolve(@ZeroFunction, acThetaStart, options);
        k_ac = GetKacFromAngles(acTheta, acFreq);
        acInvWavelen = get_angles_from_vector(k_ac) / 2 / pi;
        k_d_implied = k_ac + k_i;
        [k_d_implied_norm, dTheta, dPhi] = get_angles_from_vector(k_d_implied);
        nOrd = teo2.find_n_op( dTheta );
        k_d_norm = 2 * pi * nOrd / opWavelenVac;
        
%         % For checking this is working when results are unexpected
%         figure();
%         if isscalar(iTheta)
%             PlotWavevectorTriangle();
%         else
%             PlotErrors();
%         end
        
        function f = ZeroFunction(acTheta)
            k_ac_trial = GetKacFromAngles(acTheta, acFreq);
            k_d_implied_trial = k_ac_trial + k_i;
            [k_d_implied_norm, k_d_theta, ~] = get_angles_from_vector(k_d_implied_trial);
            nOrd_trial = teo2.find_n_op( k_d_theta );
            k_d_trial_norm = 2 * pi * nOrd_trial / opWavelenVac;
            f = k_d_implied_norm - k_d_trial_norm;
        end
        
        function PlotWavevectorTriangle()
            k_d = get_vector_from_angles(k_d_norm, dTheta, dPhi);
            nat_k_d = RotateToNaturalAxes(k_d);
            nat_k_i = RotateToNaturalAxes(k_i);
            nat_k_ac = RotateToNaturalAxes(k_ac);
            line([0 nat_k_i(1)],[0 nat_k_i(2)], [0 nat_k_i(3)], 'color', 'red'); % incident
            line([0 nat_k_d(1)],[0 nat_k_d(2)], [0 nat_k_d(3)], 'color', 'blue'); % diffracted
            line([nat_k_i(1) nat_k_i(1)+nat_k_ac(1)], [nat_k_i(2) nat_k_i(2)+nat_k_ac(2)], ...
                [nat_k_i(3) nat_k_i(3)+nat_k_ac(3)], 'color', 'green');
            grid on;
            xlabel('[110]');
            ylabel('[~10]');
            zlabel('[001]');
        end
        
        function PlotErrors()
            errorSize = k_d_implied_norm - k_d_norm;
            plot(abs(errorSize));
            xlabel('index');
            ylabel('error');
        end
    end
end

function k_ac = GetKacFromAngles(acTheta, acFreq)
    acPhi = pi/4; % Only diffraction along z axis.
    k_ac_norm = 2 * pi * acFreq ./ teo2.find_v_ac_min(acTheta, acPhi + 0*acTheta);
    k_ac = get_vector_from_angles(k_ac_norm, acTheta, acPhi);
end

function rotatedVector = RotateToNaturalAxes( vector )
    rotatedVector(1,:) = ( vector(1,:) + vector(2,:) )/2;
    rotatedVector(2,:) = ( -vector(1,:) + vector(2,:) )/2;
    rotatedVector(3,:) = vector(3,:);
end

