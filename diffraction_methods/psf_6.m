function psf_6(time, w2f, w2s, w3, w4, ws)
    a1 = [0, w2f-w2s, w3, w4, 0];
    a2 = [0, w2f-w2s, w3, w4, 0];    
    a3 = [0, w2f-w2s, w3, w4, 0];
    a4 = [0, w2f+w2s, w3, w4, 0];    
    a5 = [0, w2f+w2s, w3, w4, 0];    
    a6 = [0, w2f+w2s, w3, w4, 0];    
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {a1, a2, a3, a4, a5, a6};
    
    a = aol_fft();
    a.z_list = a.z_list + 0e-6;
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    a.xy_xz_plot_psf(propagated_wave, x, y, 0.5);
end
