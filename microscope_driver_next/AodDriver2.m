classdef AodDriver2
    properties
        params    
    end
    
    methods
        function d = AodDriver2(params)
            d.params = params;
        end
        
        function [start, stop] = convert_normalised_to_cartesian(d, scan)           
            start_z = d.remove_z_normalisation(scan.start_norm, scan);
            stop_z  = d.remove_z_normalisation(scan.stop_norm, scan);
            start_xy = d.remove_xy_normalisation(scan.start_norm, start_z, scan); 
            stop_xy  = d.remove_xy_normalisation(scan.stop_norm, stop_z, scan);
            
            start = [start_xy; start_z];
            stop  = [stop_xy; stop_z];
        end
        
        function z = remove_z_normalisation(d, normalised, scan) 
            z_norm = normalised(3,:);
            z_norm_non_zero = (z_norm == 0) * 1e-6 + z_norm;
            z = d.params.aod_aperture(end) ./ (4 * scan.acceptance_angle .* z_norm_non_zero);                  
        end
        
        function xy = remove_xy_normalisation(d, normalised, z, scan)
            xy_norm = normalised(1:2,:);
            xy = zeros(size(xy_norm));
            for m = 1:2
                xy(m,:) = 2 * scan.acceptance_angle .* z .* xy_norm(m,:); % 2 AODs -> factor of 2
            end
        end
        
        function drive_coeffs_transducer = shift_coeffs_for_transducers(d, drive_coeffs_aod_centres)
            drive_coeffs_transducer = zeros(size(drive_coeffs_aod_centres));
            time_shifts = d.params.transducer_centre_dist ./ d.params.aod_ac_vel;
            for n = 1:d.params.num_aods
                drive_coeffs_transducer(:,:,n) = shift_poly(drive_coeffs_aod_centres(:,:,n), time_shifts(n));
            end
        end
        
        function drive_coeffs_aod_centres = calc_drive_coefficients(d, start, stop, scan)
            % structure is drive_coeffs(4, num_drives, num_aods) 
            
            if d.params.num_aods == 4
                drive_coeffs_aod_centres = calc_drive_4(d.params, start, stop, scan);
            else
                error('not implemented...')
            end
        end
    end
end