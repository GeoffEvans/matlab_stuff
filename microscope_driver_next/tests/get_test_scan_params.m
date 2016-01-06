function a = get_test_scan_params()

a = ScanParams();
a.start_norm = [0, 0, 1];
a.stop_norm = [0, 0, 1];
a.ramp_time = 2e-6;

a.aod_opt_freq = [40e6, 40e6, 40e6, 40e6];
a.acceptance_angle = 4.35e-3;
a.op_wavelength = 800e-9;
a.waves_of_4th = [0,0];

a.zoom_factor = 1;
a.voxel_density_1D = 100;

end

