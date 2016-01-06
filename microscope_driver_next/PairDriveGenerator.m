classdef PairDriveGenerator
    properties
        lambda
        V
        L
        F0
        z_start
        z_stop
        x_start
        x_stop
        time_half
        D
    end
    methods
        function p = PairDriveGenerator(op_wavelength, ac_vel, pair_separation,...
                first_aod_opt_freq, x_start, x_stop, z_start, z_stop, ramp_time, D)
            p.lambda = op_wavelength;
            p.V = ac_vel;
            p.L = pair_separation;
            p.F0 = first_aod_opt_freq;
            p.z_start = z_start;
            p.z_stop = z_stop;
            p.x_start = x_start;
            p.x_stop = x_stop;
            p.time_half = ramp_time/2;
            p.D = D;
        end
        
        function drive_coeffs = get_coefficients_at_centre(p, centres_to_ref)
            [drive1, drive2] = p.get_coefficients_at_ref();
            
            drive_coeffs = zeros(4, size(drive1,2), 2);
            drive_coeffs(:,:,1) = shift_poly(drive1, centres_to_ref(1) ./ p.V);
            drive_coeffs(:,:,2) = shift_poly(drive2, -centres_to_ref(2) ./ p.V);
        end
        
        function [drive1, drive2] = get_coefficients_at_ref(p)
            s = WavefrontPair(p.lambda, p.V, p.L, p.F0);
            [w1, w2] = s.calc_wavefronts_at_ref(p.x_start, p.x_stop, p.z_start, p.z_stop, p.time_half, p.D);
            
            [drive1, drive2] = p.convert_to_drive_coeffs(w1, w2);
            
            shift = p.lambda .* p.F0 ./ p.V .* p.L; % traversed by base ray
            drive1 = shift_poly(drive1, shift ./ p.V);
        end
        
        function [drive1, drive2] = convert_to_drive_coeffs(p, w1, w2)
            % assumes -1 mode
            a1 = (p.V .^ 1 ./ p.lambda) ./ 1 .* w1(1,:);
            b1 = (p.V .^ 2 ./ p.lambda) ./ 1 .* w1(2,:) .* -1;
            c1 = (p.V .^ 3 ./ p.lambda) ./ 2 .* w1(3,:);
            d1 = (p.V .^ 4 ./ p.lambda) ./ 6 .* w1(4,:) .* -1;
            
            a2 = (p.V .^ 1 ./ p.lambda) ./ 1 .* w2(1,:) .* -1;
            b2 = (p.V .^ 2 ./ p.lambda) ./ 1 .* w2(2,:) .* -1;
            c2 = (p.V .^ 3 ./ p.lambda) ./ 2 .* w2(3,:) .* -1;
            d2 = (p.V .^ 4 ./ p.lambda) ./ 6 .* w2(4,:) .* -1;
            
            drive1 = [a1; b1; c1; d1];
            drive2 = [a2; b2; c2; d2];
        end
    end
end