function psf_4()
tic
    time = 0e-6;
    x1 = [0, 1, 0, 0, 0];
    y1 = [0, 1, 0, 0, 0];    
    x2 = [0, 1, 0, 0, 0];
    y2 = [0, 1, 0, 0, 0];    
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = 0;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    a = aol_fft();
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    a.xy_xz_plot_psf(propagated_wave, x, y, 0.5, 0.5);
toc
end

