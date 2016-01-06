function [ a, b, c, d, ticks_per_ramp ] = calc_drive_for_lv_wrapper(scan_var, system_var, xyz_norm)

aol_params = AolParams(system_var);
scan_params = ScanParams(scan_var, xyz_norm);

drive_coeffs_transducers = calc_drive_for_lv(aol_params, scan_params);
% num_orders, num_drives, num_aod)

[a, b, c, d, ticks_per_ramp] = aol_params.get_params_for_lv(drive_coeffs_transducers, scan_params.ramp_time);

end