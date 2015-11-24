function psf_4_(aol, time, a, z, v, w3, w4, ws, wf, slice)
    z_x = z + 0.05; % TODO check whether this is the right partitioning of TeO2
    L_x = 0.1;
    z_y = z;
    L_y = 0.1;    
    b = - aol.V.^2 .* aol.k ./ (2*pi) .* [...
        (1 + v(1)) ./ (L_x .* (1 + v(1)) + 2 * z_x);...
        (1 + v(2)) ./ (L_y .* (1 + v(2)) + 2 * z_y);... 
        (1 - v(1)) ./ (2 * z_x);...
        (1 - v(2)) ./ (2 * z_y)];
    
    w1 = a * (aol.half_width ./ aol.V)^1 / 1;
    w2 = b * (aol.half_width ./ aol.V)^2 / 2;

    x1 = [w1(1), w2(1), w3, w4, 0]; % 36 waves of w2f gives 1 m focal length or 64 um after obj. 
    y1 = [w1(2), w2(2), w3, w4, 0];    
    x2 = [w1(3), w2(3), w3, w4, 0];
    y2 = [w1(4), w2(4), w3, w4, 0];   
    
    waves = struct();
    waves.focus = wf;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    [propagated_wave, x, y] = aol.calculate_psf_through_aol(waves, time);
    aol.xy_xz_plot_psf(propagated_wave, x, y, slice);
    %get_psf_dimensions(aol, propagated_wave, x, y, aol.z_list)
end
