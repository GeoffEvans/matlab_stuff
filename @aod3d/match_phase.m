function [ dTheta, dPhi, k_ac_offset_norm ] = match_phase( iTheta, iPhi, acFreq, acTheta )
% Finds the diffracted and acoustic beam angles for arrays of incident
% angles and acoustic frequencies.

% Error checking: all inputs same length arrays
if sum(0*iTheta ~= 0*iPhi & 0*iTheta ~= 0*acFreq & 0*iTheta ~= 0*acTheta) > 0
    error('Check phase match params');
    %   only hoz arrays
    % 0 < Theta < pi
end

acPhi = pi/4;
opWavelenVac = aod3d.opWavelenVac;
k_i = GetKi(iTheta,iPhi);
k_ac = GetKac(acTheta, acFreq);
[dTheta, dPhi, k_ac_offset_norm] = SolveVectorTriangle(acTheta);

    function [dTheta, dPhi, k_ac_offset_norm] = SolveVectorTriangle(acTheta)      
        initialOffset = 0*acTheta;
        k_ac_offset_direction = get_vector_from_angles(1,acTheta + pi/2, acPhi+acTheta*0);
                
        options = optimset('Algorithm','trust-region-reflective','Display','off','JacobPattern',speye(length(acTheta)));
        k_ac_offset_norm = fsolve(@ZeroFunction, initialOffset, options) * 1e5;
        
        k_ac_offset = ScaleVector(k_ac_offset_direction, k_ac_offset_norm);
        k_d = k_ac + k_i + k_ac_offset;
        [~, dTheta, dPhi] = get_angles_from_vector(k_d);
        
        %MakeErrorPlot(isscalar(acTheta));
        
        function f = ZeroFunction(offsetNorm)
            k_ac_offset_trial = ScaleVector(k_ac_offset_direction, offsetNorm*1e5);
            k_d_trial = k_ac + k_i + k_ac_offset_trial;
            [k_d_norm, k_d_theta, ~] = get_angles_from_vector(k_d_trial);
            nOrd_trial = teo2.find_n_op( k_d_theta );
            k_d_required_norm = 2 * pi * nOrd_trial / opWavelenVac;
            f = k_d_norm - k_d_required_norm;
        end
                
        function MakeErrorPlot(singletons)
            figure();
            nOrd = teo2.find_n_op( dTheta );
            k_d_required_norm = 2 * pi * nOrd / opWavelenVac;
            if singletons
                PlotWavevectorTriangle();
            else
                PlotErrors();
            end
           	function PlotWavevectorTriangle()
                k_d_fromAngles = get_vector_from_angles(k_d_required_norm, dTheta, dPhi);
                nat_k_d = RotateToNaturalAxes(k_d_fromAngles);
                nat_k_i = RotateToNaturalAxes(k_i);
                nat_k_ac = RotateToNaturalAxes(k_ac);
                line([0 nat_k_i(1)],[0 nat_k_i(2)], [0 nat_k_i(3)], 'color', 'red'); % incident
                line([0 nat_k_d(1)],[0 nat_k_d(2)], [0 nat_k_d(3)], 'color', 'blue'); % diffracted
                line([nat_k_i(1) nat_k_i(1)+nat_k_ac(1)], [nat_k_i(2) nat_k_i(2)+nat_k_ac(2)], ...
                    [nat_k_i(3) nat_k_i(3)+nat_k_ac(3)], 'color', 'green');
                nat_k_ac = RotateToNaturalAxes(k_ac + k_ac_offset);
                line([nat_k_i(1) nat_k_i(1)+nat_k_ac(1)], [nat_k_i(2) nat_k_i(2)+nat_k_ac(2)], ...
                    [nat_k_i(3) nat_k_i(3)+nat_k_ac(3)], 'color', 'black');
                grid on;
                xlabel('[110]');
                ylabel('[~10]');
                zlabel('[001]');
                function rotatedVector = RotateToNaturalAxes( vector )
                    rotatedVector(1,:) = ( vector(1,:) + vector(2,:) )/2;
                    rotatedVector(2,:) = ( -vector(1,:) + vector(2,:) )/2;
                    rotatedVector(3,:) = vector(3,:);
                end
            end
            function PlotErrors()
                [k_d_norm, ~, ~] = get_angles_from_vector(k_d);
                errorSize = k_d_norm - k_d_required_norm;
                plot(abs(errorSize));
                xlabel('index');
                ylabel('error');
            end
        end
    end
    function k_ac = GetKac(acTheta, acFreq)
        k_ac_norm = 2 * pi * acFreq ./ teo2.find_v_ac_min(acTheta, acPhi+acTheta*0);
        k_ac = get_vector_from_angles(k_ac_norm, acTheta, acPhi+acTheta*0);
    end
    function k_i = GetKi(iTheta,iPhi)
        [ ~, nExt, ~, ~ ] = teo2.find_n_op( iTheta );
        k_i_norm = 2 * pi * nExt / opWavelenVac;
        k_i = get_vector_from_angles(k_i_norm, iTheta, iPhi);
    end
    function v = ScaleVector(vIn, factor)
        v(1,:) = vIn(1,:) .* factor;
        v(2,:) = vIn(2,:) .* factor;
        v(3,:) = vIn(3,:) .* factor;
    end
end


