function drive_coeffs_transducers = calc_drive_for_lv(aol_params, scan_params)

driver = AodDriver2(aol_params);
[start, stop] = driver.convert_normalised_to_cartesian(scan_params);
drive_coeffs_aod_centres = driver.calc_drive_coefficients(start, stop, scan_params);
drive_coeffs_transducers = driver.shift_coeffs_for_transducers(drive_coeffs_aod_centres);

end
