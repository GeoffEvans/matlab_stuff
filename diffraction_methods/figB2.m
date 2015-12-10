aol = aol_fft();
aol.adjustment = 200;
aol.number_of_samples = 2^8 - 1;
aol.z_list = linspace(-40,40,499)*1e-6;
aol.cm = cm;
aol.spacing = 0;

w_list = linspace(-5,5,41);
z_fwhm0_list = [];
z_fwhm4_list = [];
z_fwhm6_list = [];
x_fwhm0_list = [];
x_fwhm4_list = [];
x_fwhm6_list = [];

for ws = w_list
   [x_fwhm0, z_fwhm0] = psf_4(aol, 0, [0,0,0,0], 0, 0, 0, 0, 0, ws, 0, 0.5);
   [x_fwhm4, z_fwhm4] = psf_4(aol, 0, [0,0,0,0], 0, 0, 0, -ws/2, 0, ws, 0, 0.5);
   [x_fwhm6, z_fwhm6] = psf_6(aol, 0, 0, 0, 0, 0, -ws/4*1.8, 0, ws, 0, 0.5);
   x_fwhm0_list = [x_fwhm0_list, x_fwhm0];
   x_fwhm4_list = [x_fwhm4_list, x_fwhm4];
   x_fwhm6_list = [x_fwhm6_list, x_fwhm6];
   z_fwhm0_list = [z_fwhm0_list, z_fwhm0];
   z_fwhm4_list = [z_fwhm4_list, z_fwhm4];
   z_fwhm6_list = [z_fwhm6_list, z_fwhm6];
end
hold on
qq = linspace(-4,4);
plot(qq, spline(w_list, z_fwhm0_list, qq), 'r')
plot(qq, spline(w_list, z_fwhm4_list, qq), 'b')
plot(qq, spline(w_list, z_fwhm6_list, qq), 'k')
plot(qq, spline(w_list, x_fwhm0_list, qq), 'r--')
plot(qq, spline(w_list, x_fwhm4_list, qq), 'b--')
plot(qq, spline(w_list, x_fwhm6_list, qq), 'k--')

