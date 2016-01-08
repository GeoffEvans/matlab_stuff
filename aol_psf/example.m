aol = aol_fft();
aol.z_list = linspace(115, 155, 199) * 1e-6;
%psf_4_waves(aol, 0, 0, 0, 0, 0, 0, 0, 0, 0)
psf_4_drives(aol, 0, [0,0,0,0], -1, [0,0], 0, [0,0,0,0], 0, 0, 0)
%psf_6_waves(aol, 0, 0, 0, 0, 0, 0, 0, 0, 0)
%psf_6_drives(aol, 0, [0,0,0,0,0,0], 1e9, 0, 0, [0,0,0,0,0,0], 0, 0, 0)