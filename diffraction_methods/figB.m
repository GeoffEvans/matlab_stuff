ws = -10;
time = (-1:1:1)*0.25e-6;
aol = aol_fft();
aol.z_list = linspace(-50,50,199)*1e-6;
aol.pwr = 1;
aol.cm = cm; % jet;

psf_6(aol, time, 0, 0, 0, 0, 0, 0.5) % no aberration
psf_6(aol, time, 0, 0, 0, 0, ws, 0.39) % spherical aberration
psf_4(aol, time, 0, 0, 0, 5, ws, 0.47) % correction using 4-AOD
psf_6(aol, time, 0, 0, 0, 4.7, ws, 0.5) % correction using 6-AOD
