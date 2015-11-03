ws = -10;
t = (-1:1:1)*5e-6;
psf_6(t, 0, 0, 0, 0, 0) % no aberration
psf_6(t, 0, 0, 0, 0, ws) % spherical aberration
psf_4(t, 0, 0, 0, 7, ws) % instantaneous correction using 4-AOD
psf_6(t, 0, 0, 0, 5, ws) % instantaneous correction using 6-AOD
