function psf_4(aol, time, a, w2f, w2s, w3, w4, ws, wf, slice)
    w1 = a * (aol.half_width ./ aol.V)^1 / 1;    

    x1 = [w1(1), w2f-w2s, w3, w4, 0]; % 36 waves of w2f gives 1 m focal length or 64 um after obj. 
    y1 = [w1(2), w2f    , w3, w4, 0];    
    x2 = [w1(3), w2f+w2s, w3, w4, 0];
    y2 = [w1(4), w2f    , w3, w4, 0];   
    
    waves = struct();
    waves.focus = wf;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    [propagated_wave, x, y] = aol.calculate_psf_through_aol(waves, time);
    aol.xy_xz_plot_psf(propagated_wave, x, y, slice);
    get_psf_dimensions(aol, propagated_wave, x, y, aol.z_list)
end
