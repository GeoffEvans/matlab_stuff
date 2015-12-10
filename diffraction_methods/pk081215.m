aol = aol_fft();
aol.z_list = linspace(-150, -110, 199) * 1e-6;

w4 = linspace(0,0,6);
z = 1; % focal length of AOL in metres
psf_6_(aol, [0 4e-6], [0,0,0,0,0,0], z, -1, 0, w4, 0, 0, 0, 0.5)
%psf_6_(aol, time, a, z, v, w3, w4, w5, ws, wf, slice)
% a is centre freq