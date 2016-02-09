aol = aol_fft();
aol.adjustment = 380;
aol.number_of_samples = 2^9 - 1;
aol.k = 2*pi/920e-9;
aol.spacing = 0;

[x_fwhm6, z_fwhm6, max_fl6, total_fl6] = plot_na(aol, 0.6);
[x_fwhm7, z_fwhm7, max_fl7, total_fl7] = plot_na(aol, 0.7);
[x_fwhm8, z_fwhm8, max_fl8, total_fl8] = plot_na(aol, 0.8);

% aol.k = 2*pi/800.1e-9;
% 
% [x_fwhm6c, z_fwhm6c] = plot_na(aol, 0.6);
% [x_fwhm8c, z_fwhm8c] = plot_na(aol, 0.8);