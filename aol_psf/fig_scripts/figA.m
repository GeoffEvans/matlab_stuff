time = (-10:1:10) * 1e-6;
aol = aol_fft();
aol.z_list = linspace(-20,20,299)*1e-6;
aol.cm = cm; % jet;

psf_4(aol, time, 0, 4, 0, 0, 0, 0.5); % x
psf_4(aol, time, 0, 4, 0.8, 0, 0, 0.5); % xz
psf_4(aol, time, 1.3, 4, 0, -0.5, 0, 0.5); % curved