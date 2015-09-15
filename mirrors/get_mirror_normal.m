function normal = get_mirror_normal( k1, k2 )
    normal = k1 - k2;
    normal = normal / norm(normal);
    if normal(1) < 0
        normal = -normal;
    end
end

