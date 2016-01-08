aol = aol_fft();
aol.adjustment = 200;
aol.number_of_samples = 2^9 - 1;

z_list = [-250:50:-50 1e-9 50:50:250];
z_fwhm0_list = [];
x_fwhm0_list = [];

for z = z_list
   aol.z_list = linspace(z-25,z+25,199) * 1e-6;
   z_aol = -130 / z;
   res = psf_4_drives(aol, 0, [0,0,0,0], z_aol, [0,0], 0, [0,0,0,0], 0, 0, 0)
   x_fwhm0_list = [x_fwhm0_list, res(3)];
   z_fwhm0_list = [z_fwhm0_list, res(4)];
end
figure()
hold on
plot(z_list, z_fwhm0_list, 'r')
plot(z_list, x_fwhm0_list, 'b--')


