ws = -10;
time = 0;%(-1:1:1)*0.25e-6;
aol = aol_fft();
aol.adjustment = 2000;
aol.number_of_samples = 2^9 - 1;
aol.z_list = linspace(-40,10,399)*1e-6;
aol.pwr = 4;
aol.cm = cm; % jet;

psf_6(aol, time, 0, 0, 0, 0, 0, 0, -2, 0, 0.5) % no aberration
%psf_6(aol, time, 0, 0, 0, 0, 0, 0, ws, 0, 0.39) % spherical aberration
%psf_4(aol, time, [0,0,0,0], 0, 0, 0, 5, 0, ws, 0, 0.47) % correction using 4-AOD
%psf_6(aol, time, 0, 0, 0, 0, 4.7, 0, ws, 0, 0.5) % correction using 6-AOD
