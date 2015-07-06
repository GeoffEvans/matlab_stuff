function [ r_res, z_res ] = get_psf_dimensions( field_3d, x, y, z)
    r = sqrt(x.^2 + y.^2);
    max_intensity = max(field_3d(:).^2);
    [row, col, depth] = find(field_3d == max_intensity);
    display(2 * [row/size(field_3d,1), col/size(field_3d,2), depth/size(field_3d,3)]) % should print 1 1 1
    
    half_intensity_or_less_r = max(field_3d, 3) < max_intensity;
    r_half = r(half_intensity_or_less_r);
    r_res = min(r_half(:));
    
    half_intensity_or_less_z = max(max(field_3d)) < max_intensity;
    z_half = z(half_intensity_or_less_z);
    z_res = min(z_half(:));
end

