function a = get_test_aol_params()

a = AolParams();
a.num_aods = 4;
a.aol_config = AolConfig.Orig_4;
a.system_clock_freq = 200e6;
a.data_time_interval = 50e-9;

a.aod_ac_vel = repmat(613, 1, 4);
a.aod_aperture = repmat(20e-3, 1, 4);
a.transducer_centre_dist = repmat(10e-3, 1, 4);

a.aod_mode = repmat(-1, 1, 4);
a.aod_xy_offsets = -[[0;0], [0.0025;0], [0.005;0.0025], [0.005;0.005]];
a.aod1_centre_to_ref = -[0.0052, 0.0052];
a.aod_reduced_z = [-15e-2, -10e-2, -5e-2, 0];

end

