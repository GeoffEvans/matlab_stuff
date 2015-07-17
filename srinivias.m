function [corners, starts, stops] = srinivias( point1, point2, point3 )

u13 = point1 - point3;
u32 = point3 - point2;
u21 = point2 - point1;

scaling = -u13(3) ./ u32(3);
xy_displacement = u13 + scaling .* u32;
v_xy = xy_displacement ./ norm(xy_displacement);
xz_displacement = u13 - dot(u13, v_xy) .* v_xy;
v_xz = xz_displacement ./ norm(xz_displacement);

v_normal = cross(v_xy, v_xz);
normal_shift = repmat(v_normal, 1, 4) * dot(point1, v_normal);

vectors = {point1, point2, point3};

pad = 4;
min_xy = min(cellfun(@(v) dot(v_xy, v), vectors)) - pad / norm(v_xy(1:2));
max_xy = max(cellfun(@(v) dot(v_xy, v), vectors)) + pad / norm(v_xy(1:2));
min_xz = min(cellfun(@(v) dot(v_xz, v), vectors)) - pad / norm(v_xz(1:2));
max_xz = max(cellfun(@(v) dot(v_xz, v), vectors)) + pad / norm(v_xz(1:2));

corners = [ min_xy * v_xy + min_xz * v_xz,...
            min_xy * v_xy + max_xz * v_xz,...
            max_xy * v_xy + max_xz * v_xz,...
            max_xy * v_xy + min_xz * v_xz ] + normal_shift;

number_of_lines = round((max_xz - min_xz) .* norm(v_xz(1:2)));
frac = linspace(0, 1, number_of_lines);

starts = corners(:,1) * frac + corners(:,2) * (1 - frac);
stops = corners(:,4) * frac + corners(:,3) * (1 - frac);
        
end

