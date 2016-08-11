V = 613;
aperture = 16e-3;
scan_time = 20e-6;
spacing = 0.0;
focus = -0.8167;
axial_scan = -39.5;

axial_range = -2;
num_ray = 11;
num_times = 1;
[rays_t, rays_x] = meshgrid(linspace(-scan_time/2, scan_time/2, num_times), linspace(-aperture/2, aperture/2, num_ray));

B = 1 + spacing/2/focus;
w3 = axial_scan / focus / B^3
w4 = - 3 * spacing * B * w3^2
w5 = 15 * spacing^2 / B * w3^3
phase_aod_first = @(t, x) (x - V*t).^2 / 2 / (2 * focus + spacing) + axial_scan / focus * (x - V*t).^3 / 6; % + w4 * (x - V*t).^4 / 24 + w5 * (x - V*t).^5 / 120;
phase_aod_second = @(t, x) (x + V*t).^2 / 2 / (2 * focus) - (x + V*t).^3 / 6 * axial_scan / focus; 

D = 1e-9;
deflect = @(t, x, phase_aod) -(phase_aod(t, x+D) - phase_aod(t, x-D)) ./ (2*D);
propagate = @(k, L) k * L;

rays_k = deflect(rays_t, rays_x, phase_aod_first);
rays_x = rays_x + propagate(rays_k, spacing);
rays_k = rays_k + deflect(rays_t, rays_x, phase_aod_second);

X = [rays_x(:), rays_x(:) + propagate(rays_k(:), axial_range)];
Y = [zeros(num_ray * num_times, 1), zeros(num_ray * num_times, 1) + axial_range];
figure()
plot(X', Y');


