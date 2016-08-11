function wavefronts_test()

    V = 613;
    aperture = 16e-3;
    scan_time = 2e-6;
    spacing = 0.1;
    focus = 2;
    axial_scan = 900;
    axial_range = 4;
    deflection1 = 0e-3;
    deflection2 = 0e-3;

    B = 1 + spacing/2/focus;
    w3 = axial_scan / focus / B^3;
    phase_aod_first = @(t, x) (x - V*t)*deflection1 + (x - V*t).^2 / 2 / (2 * focus + spacing) + w3 * (x - V*t).^3 / 6 - 3 * spacing * B * w3^2 * (x - V*t).^4 / 24 + 15 * spacing^2 / B * w3^3 * (x - V*t).^5 / 120;
    phase_aod_second = @(t, x) (x + V*t)*deflection2 + (x + V*t).^2 / 2 / (2 * focus) - (x + V*t).^3 / 6 * axial_scan / focus; 
    
    figure()
    hold on;
    
    plot_rays(axial_range, scan_time, aperture, spacing, phase_aod_first, phase_aod_second)
    plot_phase(axial_range, 920e-9, scan_time, aperture, spacing, phase_aod_first, phase_aod_second)
end
    
function plot_rays(axial_range, scan_time, aperture, spacing, phase_aod_first, phase_aod_second)
    num_ray = 11;
    num_times = 1;
    [rays_t, rays_x] = meshgrid(linspace(-scan_time/2, scan_time/2, num_times), linspace(-aperture/2, aperture/2, num_ray));

    D = 1e-9;
    deflect = @(t, x, phase_aod) -(phase_aod(t, x+D) - phase_aod(t, x-D)) ./ (2*D);
    propagate = @(k, L) k * L;

    rays_k = deflect(rays_t, rays_x, phase_aod_first);
    rays_x = rays_x + propagate(rays_k, spacing);
    rays_k = rays_k + deflect(rays_t, rays_x, phase_aod_second);

    X = [rays_x(:), rays_x(:) + propagate(rays_k(:), axial_range)];
    Y = [zeros(num_ray * num_times, 1), zeros(num_ray * num_times, 1) + axial_range];
    plot(X', Y');
end

function plot_phase(axial_range, wavelength, scan_time, aperture, spacing, phase_aod_first, phase_aod_second)
    wavevector = 2 * pi / wavelength;
    width = aperture*11;
    fft_number_of_samples = 2^14 - 1;
    z_list = 2.88 + linspace(-0.05, 0.05, 111);
    x_list = linspace(-width/2, width/2, fft_number_of_samples)';
    
    wave = abs(x_list) < aperture/2;
    wave = wave ./ exp(1i * wavevector * phase_aod_first(scan_time/2, x_list)); % first AOD phase
    wave = propagate_wave(wave, wavevector, spacing, width, fft_number_of_samples, 1); % propagate to second
    wave = wave ./ exp(1i * wavevector * phase_aod_second(scan_time/2, x_list)); % combine with second AOD phase
    wave_range = propagate_wave(wave, wavevector, z_list, width, fft_number_of_samples, 1); % propagate through focal region
    
    [X, Z] = meshgrid(x_list, z_list);
    pcolor(X, Z, angle(wave_range)')
    
    calculate_psf(wave, width, fft_number_of_samples, wavevector);
end

function propagated_wave = propagate_wave(wave, wavevector, distances, space_width, fft_number_of_samples, ref_ind)
    count = [0:floor((fft_number_of_samples - 1)/2) floor(-(fft_number_of_samples - 1)/2):-1]';

    k_x = count * 2 * pi / space_width;
    k_z = sqrt(wavevector.^2 * ref_ind.^2 - k_x.^2);
    mask = ~(abs(imag(k_z)) > 0);

    ft_wave = fft(ifftshift(wave)); % recentre wave with fftshift

    propagated_wave = zeros([size(k_z, 1), numel(distances)]);
    for n = 1:numel(distances) % to handle multiple distances in one hit
        propagated_wave(:,n) = fftshift(ifft((ft_wave .* exp(1i * mask .* k_z .* distances(n)) .* mask))); % shift back 
    end
    propagated_wave = squeeze(propagated_wave); % in case only one distance required
end

function calculate_psf(wave_leaving_aol, width, fft_number_of_samples, wavevector)
    tube_lens_focal_length = 160e-3;  
    aol_to_objective_scaling = 0.8;
    objective_magnification = 20;
    
    integer_list = floor(-fft_number_of_samples/2+0.5):floor(fft_number_of_samples/2-0.5);
    k_x_discrete = integer_list * 2 * pi / (width * aol_to_objective_scaling);
    plane_wave_angle_air = k_x_discrete / wavevector;
    x = plane_wave_angle_air * tube_lens_focal_length / objective_magnification;
    width_fp = (max(x(:)) - min(x(:)));
    
    focal_plane_wave = fftshift(fft2(ifftshift(wave_leaving_aol))); % transform to focal plane    
    
    z_range = linspace(-1,1,101)*20e-6;
    focal_region_wave = propagate_wave(focal_plane_wave, wavevector, z_range, width_fp, fft_number_of_samples, 4/3);
    
    [X,Z] = meshgrid(x, z_range);
    
    figure()
    pcolor(X, Z, abs(focal_region_wave)');
    axis equal
end