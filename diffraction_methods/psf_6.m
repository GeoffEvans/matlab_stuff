function psf_6(aol, time, w2f, w2s, w3, w4, ws, wf, slice)
    a1 = [0, w2f-w2s, w3, w4, 0];
    a2 = [0, w2f    , w3, w4, 0];    
    a3 = [0, w2f    , w3, w4, 0];
    a4 = [0, w2f+w2s, w3, w4, 0];    
    a5 = [0, w2f    , w3, w4, 0];    
    a6 = [0, w2f    , w3, w4, 0];    
    
    waves = struct();
    waves.focus = wf;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {a1, a2, a3, a4, a5, a6};
    
    [propagated_wave, x, y] = aol.calculate_psf_through_aol(waves, time);
    aol.xy_xz_plot_psf(propagated_wave, x, y, slice);
    get_psf_dimensions(aol, propagated_wave, x, y, aol.z_list)
end
