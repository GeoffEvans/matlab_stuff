function [x_fwhm, z_fwhm] = psf_6(aol, time, w1, w2f, w2s, w3, w4, w5, ws, wf, slice)
    a1 = [w1, w2f-w2s, w3, w4, w5];
    a2 = [w1, w2f    , w3, w4, w5];    
    a3 = [w1, w2f    , w3, w4, w5];
    a4 = [w1, w2f+w2s, w3, w4, w5];    
    a5 = [w1, w2f    , w3, w4, w5];    
    a6 = [w1, w2f    , w3, w4, w5];    
    
    waves = struct();
    waves.focus = wf;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {a1, a2, a3, a4, a5, a6};
    
    [propagated_wave, x, y] = aol.calculate_psf_through_aol(waves, time);
    aol.xy_xz_plot_psf(propagated_wave, x, y, slice);
    res = get_psf_dimensions(aol, propagated_wave, x, y, aol.z_list);
    x_fwhm = res(3);
    z_fwhm = res(4);
end
