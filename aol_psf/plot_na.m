function [x_fwhm, z_fwhm, max_fl, total_fl] = plot_na(aol, NA)

aol.half_width = NA * 1e-2;
aol.beam_sigma = NA/2 * 1e-2; 

z_res = 199; 

z_list = [1e-9 50:50:250];
x_list = 0:50:250;
[x_mesh, z_mesh] = meshgrid(x_list, z_list);
points = size(x_mesh);
z_fwhm = zeros(points);
x_fwhm = zeros(points);
max_fl = zeros(points);
total_fl = zeros(points);

for m = 1:points(1)
    for n = 1:points(2)
        z = z_mesh(m,n);
        x = x_mesh(m,n);
        fprintf('%f %f\n', [x z])
        aol.z_list = linspace(z-25,z+25,z_res) * 1e-6;
        z_aol = -130 / z;
        x_aol = z_aol * x * 1e-4; % assumes 0.8 mag relay
        res = psf_4_linear(aol, 0, [x_aol, 0], z_aol, [0,0], 0, [0,0,0,0], 0, 0, 0, false);
        x_fwhm(m,n) = res(3);
        z_fwhm(m,n) = res(4);
        max_fl(m,n) = res(5);
        total_fl(m,n) = res(6);
    end
end
x_mesh = [fliplr(-x_mesh), x_mesh; fliplr(-x_mesh), x_mesh];
z_mesh = [flipud(-z_mesh), flipud(-z_mesh); z_mesh, z_mesh];
x_fwhm = expand(x_fwhm);
z_fwhm = expand(z_fwhm);
max_fl = expand(10.^(max_fl - max(max_fl(:))));
total_fl = expand(10.^(total_fl - max(total_fl(:))));

ellipse_grid(x_mesh, z_mesh, x_fwhm, z_fwhm);
ellipse_grid(x_mesh, z_mesh, max_fl, total_fl);
end

function expanded = expand(mat)
    expanded = [rot90(mat,2), flipud(mat); fliplr(mat), mat];
end