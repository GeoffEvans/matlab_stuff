function drive_coeffs_aod_centres = calc_drive_4(aol, start, stop, scan)
    MAX_ORDER = 4; % a,b,c,d  
    drive_coeffs_aod_centres = zeros(MAX_ORDER, size(start, 2), aol.num_aods);   
    
    pairings = aol.aol_config.get_pairings();
    
    [drive_coeffs_aod_centres(:,:,pairings(1,:))] = calc_drive_pair_from_4(...
        start(1,:), stop(1,:), start(3,:), stop(3,:), aol, scan, 1);
    
    [drive_coeffs_aod_centres(:,:,pairings(2,:))] = calc_drive_pair_from_4(...
        start(2,:), stop(2,:), start(3,:), stop(3,:), aol, scan, 2);
end

function drive_pair = calc_drive_pair_from_4(start_x, stop_x, start_z, stop_z, aol, scan, x_or_y)
    pairings = aol.aol_config.get_pairings();
    pair_idx = pairings(x_or_y, :);
    
    pair_separation = aol.aod_reduced_z(pair_idx(2)) - aol.aod_reduced_z(pair_idx(1));
    z_start_eff = start_z - aol.aod_reduced_z(pair_idx(2));
    z_stop_eff  = stop_z  - aol.aod_reduced_z(pair_idx(2)); 

    centre_to_ref_offsets = repmat(aol.aod1_centre_to_ref(1:2)', 1, aol.num_aods) - aol.aod_xy_offsets;
    
    pair = PairDriveGenerator(scan.op_wavelength, aol.aod_ac_vel(pair_idx(1)),...
        pair_separation, scan.aod_opt_freq(pair_idx(1)), start_x, stop_x,...
        z_start_eff, z_stop_eff, scan.ramp_time, scan.get_D(aol, x_or_y));
    
    drive_pair = pair.get_coefficients_at_centre(centre_to_ref_offsets(x_or_y, pair_idx));
end