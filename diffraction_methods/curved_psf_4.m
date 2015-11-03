max_wave = 0;
for t = -3:0.03:3
    time = t * 1e-5;
    w1 = 0;
    w2f = 0;
    w2s = 1;
    w3 = 0.4;
    w4 = 0;
    ws = 0;
    x1 = [-w1, w2f-w2s, w3, w4, 0];
    y1 = [-w1, w2f-w2s, w3, w4, 0];    
    x2 = [w1, w2f+w2s, w3, w4, 0];
    y2 = [w1, w2f+w2s, w3, w4, 0];    
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    a = aol_fft();
    a.z_list = a.z_list + 10e-6;
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    max_wave = max(max_wave, propagated_wave);
end
a.xy_xz_plot_psf(max_wave, x, y, 0.5);    
