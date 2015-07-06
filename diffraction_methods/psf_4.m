function psf_4()
tic
    time = 0e-6;
    x1 = [0, 10, 0, 0, 0];
    y1 = [0, 10, 0, 0, 0];    
    x2 = [0, 10, 0, 0, 0];
    y2 = [0, 10, 0, 0, 0];    
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = 0;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    a = aol_fft();
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    a.xy_xz_plot_psf(propagated_wave, x, y, 0.5);
    get_psf_dimensions(propagated_wave, x, y, a.z_list)
toc
end

function res = get_psf_dimensions(field_3d, x, y, z)
    r = sqrt(x.^2 + y.^2);
    max_intensity_sqr = max(abs(field_3d(:)).^4)
    [row, col, depth] = find(field_3d == max_intensity_sqr);
    display(2 * [row/size(field_3d,1), col/size(field_3d,2), depth/size(field_3d,3)]) % should print 1 1 1
    
    half_or_more_r = max(abs(field_3d), [], 3).^4 >= max_intensity_sqr/2;
    r_res = 2 * max(r(half_or_more_r));
    
    half_or_more_z = squeeze(max(max(abs(field_3d)))).^4 >= max_intensity_sqr/2;
    z_res = max(z(half_or_more_z)) - min(z(half_or_more_z));
    
    res = [r_res, z_res] * 1e6;
end