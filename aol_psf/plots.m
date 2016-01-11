aol = aol_fft();
aol.adjustment = 100;
aol.number_of_samples = 2^9 - 1;
aol.k = 2*pi/800e-9;

z_res = 199; 
res = 66.66;

z_list = [(-200:res:-res) 1e-9 (res:res:200)];
x_list = -200:res:200;
[x_mesh, z_mesh] = meshgrid(x_list, z_list);
points = numel(x_mesh);
z_fwhm0_list = zeros(1, points);
x_fwhm0_list = zeros(1, points);

for n = 1:points
    z = z_mesh(n);
    x = x_mesh(n);
    aol.z_list = linspace(z-25,z+25,z_res) * 1e-6;
    z_aol = -130 / z;
    x_aol = z_aol * x * 1e-4; % assumes 0.8 mag relay
    res = psf_4_linear(aol, 0, [x_aol, 0], z_aol, [0,0], 0, [0,0,0,0], 0, 0, 0);
    xlabel(x)
    ylabel(z)
    x_fwhm0_list(n) = res(3);
    z_fwhm0_list(n) = res(4);
end

ellipse_grid(x_mesh(:), z_mesh(:), x_fwhm0_list', z_fwhm0_list')


