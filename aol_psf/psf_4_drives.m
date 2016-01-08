function res = psf_4_drives(aol, time, a, z, v, w3, w4, w5, ws, wf)
    z_x = z + aol.spacing;
    L_x = 2 * aol.spacing;
    z_y = z;
    L_y = 2 * aol.spacing;    
    b = - aol.V.^2 .* aol.k ./ (2*pi) .* [...
        (1 + v(1)) ./ (L_x .* (1 + v(1)) + 2 * z_x);...
        (1 + v(2)) ./ (L_y .* (1 + v(2)) + 2 * z_y);... 
        (1 - v(1)) ./ (2 * z_x);...
        (1 - v(2)) ./ (2 * z_y)];
    
    w1 = a * (aol.half_width ./ aol.V)^1 / 1;
    w2 = b * (aol.half_width ./ aol.V)^2 / 2;

    x1 = [w1(1), w2(1), w3, w4(1), w5]; % 36 waves of w2f gives 1 m focal length or 64 um after obj. 
    y1 = [w1(2), w2(2), w3, w4(2), w5];    
    x2 = [w1(3), w2(3), w3, w4(3), w5];
    y2 = [w1(4), w2(4), w3, w4(4), w5];   
    aod_drives_in_waves = {x1, y1, x2, y2};
    
    res = get_psf(aol, time, ws, wf, aod_drives_in_waves, true);
end
