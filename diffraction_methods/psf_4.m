function psf_4(time, w2f, w2s, w3, w4, ws)
    x1 = [0, w2f-w2s, w3, w4, 0];
    y1 = [0, w2f-w2s, w3, w4, 0];    
    x2 = [0, w2f-w2s, w3, w4, 0];
    y2 = [0, w2f+w2s, w3, w4, 0];   
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    a = aol_fft();
    a.z_list = a.z_list;
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    a.xy_xz_plot_psf(propagated_wave, x, y, 0.5);
end
