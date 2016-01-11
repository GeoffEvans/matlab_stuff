aol = aol_fft();
aol.adjustment = 200;
aol.k = 2*pi/800e-9;
aol.number_of_samples = 2^8 - 1;
z_res = 199; 
res = 60;

z_list = [(-240:res:-res) 1e-9 (res:res:240)];
x_list = -240:res:240;
[x_mesh, z_mesh] = meshgrid(x_list, z_list);
points = numel(x_mesh);
z_fwhm0_list = zeros(1, points);
x_fwhm0_list = zeros(1, points);

n = 1
z = z_mesh(n);
x = x_mesh(n);
aol.z_list = linspace(z-25,z+25,z_res) * 1e-6;
z_aol = -130 / z;
x_aol = z_aol * x * 1e-4; % assumes 0.8 mag relay
res = psf_4_linear(aol, 0, [x_aol, 0], z_aol, [0,0], 0, [0,0,0,0], 0, 0, 0);