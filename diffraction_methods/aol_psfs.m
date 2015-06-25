function aol_psfs()

waves_x1 = [0, 1, 0, 0, 0];
waves_y1 = [0, 1, 0, 0, 0];    
waves_x2 = [0, 1, 0, 0, 0];
waves_y2 = [0, 1, 0, 0, 0];    
time = 0e-6;

adjustment = 1e-3;
number_of_samples = 2.^10 - 1; % computational speed vs accuracy
k = 2*pi/920e-9;
z_list = linspace(-20,20,39)*1e-6;
z_plane_no = round(size(z_list,2)/2); % which z plane is plotted for the xy psf 

waves = struct('x1', waves_x1, 'y1', waves_y1, 'x2', waves_x2, 'y2', waves_y2);
do_plot = struct('input', 0, 'focal_plane', 0, 'angular_spec', 0);

calculate_and_plot(waves, time, number_of_samples, k, adjustment, z_list, z_plane_no, do_plot);
end

function calculate_and_plot(waves, time, number_of_samples, k, adjustment, z_list, z_plane_no, do_plot)
    [sampled_wave_2d, space_width] = get_sampled_wavefunction(waves, time, number_of_samples, k, do_plot.input);       
    [x, y, focal_plane_wave_2d, space_width] = get_focal_plane_wavefunction(sampled_wave_2d, adjustment, space_width, number_of_samples, k, do_plot.focal_plane);
    propagated_wave_2d = propagate_wave(focal_plane_wave_2d, z_list, space_width, k, number_of_samples, do_plot.angular_spec);
    xy_xz_plot_psf(propagated_wave_2d, x, y, z_list, z_plane_no);
end

function plot_wavefunction_2d(wave_function_2d, x, y)
    figure()
    subplot(1,2,1)
    hh = pcolor(x,y,abs(wave_function_2d));
    set(hh,'EdgeColor','none')
    axis equal
    subplot(1,2,2)
    hh = pcolor(x,y,angle(wave_function_2d));
    set(hh,'EdgeColor','none')
    axis equal
end

function [sampled_wave, space_width] = get_sampled_wavefunction(waves, t, number_of_samples, k, plot_input)     
    space_width = pi * number_of_samples / k;
    samples = linspace(-1/2, 1/2, number_of_samples) * space_width;
    [x, y] = meshgrid(samples);
    
    % take aperture to be 12mm because of scaling by relay
    scale_factors = 2*pi ./ (6e-3).^(1:5) / k;
    phase_x1 = k * (waves.x1 .* scale_factors);
    phase_y1 = k * (waves.y1 .* scale_factors);
    phase_x2 = k * (waves.x2 .* scale_factors);
    phase_y2 = k * (waves.y2 .* scale_factors);
    sph_ab_waves = -4;
    focus_waves = 3.5;    
 
    V = 613;
    sampled_wave = (sqrt(x.^2 + y.^2) < 6e-3)... 
     .* exp(1i * (phase_x1(1) .* (x-V*t) + phase_x1(2) .* (x-V*t).^2 + phase_x1(3) .* (x-V*t).^3 + phase_x1(4) .* (x-V*t).^4));
    %plot_wavefunction_2d(sampled_wave, x, y)
    sampled_wave = propagate_wave(sampled_wave, 0.05, space_width, k, number_of_samples, 0)...
     .* exp(1i * (phase_y1(1) .* (y-V*t) + phase_y1(2) .* (y-V*t).^2 + phase_y1(3) .* (y-V*t).^3 + phase_y1(4) .* (y-V*t).^4)); 
    %plot_wavefunction_2d(sampled_wave, x, y)
    sampled_wave = propagate_wave(sampled_wave, 0.05, space_width, k, number_of_samples, 0)...
     .* exp(1i * (phase_x2(1) .* (-x-V*t) + phase_x2(2) .* (-x-V*t).^2 + phase_x2(3) .* (-x-V*t).^3 + phase_x2(4) .* (-x-V*t).^4));
    %plot_wavefunction_2d(sampled_wave, x, y)
    sampled_wave = propagate_wave(sampled_wave, 0.05, space_width, k, number_of_samples, 0)...
     .* exp(1i * (phase_y2(1) .* (-y-V*t) + phase_y2(2) .* (-y-V*t).^2 + phase_y2(3) .* (-y-V*t).^3 + phase_y2(4) .* (-y-V*t).^4))...
     .* exp(1i * (x.^2 + y.^2) * focus_waves * 2*pi ./ (6e-3).^2 ./ 2)...  
     .* exp(1i * (x.^2 + y.^2).^2 * sph_ab_waves * 2*pi ./ (6e-3).^4 ./ 4);
 
    if plot_input
        plot_wavefunction_2d(sampled_wave, x, y)
    end
end

function propagated_wave = propagate_wave(wave_2d, distance, space_width, k, number_of_samples, plot_angular_spectrum)
    count = [0:floor((number_of_samples - 1)/2) floor(-(number_of_samples - 1)/2):-1];
    [k_x, k_y] = meshgrid(count * 2 * pi / space_width);
    ft_wave_2d = fft2(ifftshift(wave_2d)); % recentre wave with fftshift
       
    assert(max(k_x(:))/k > 0.99);
    k_z = sqrt(k.^2 - k_x.^2 - k_y.^2);
    mask = ~(abs(imag(k_z)) > 0);
    
    if numel(distance) == 1
        propagated_wave = fftshift(ifft2((ft_wave_2d .* exp(1i * mask .* k_z .* distance) .* mask))); % shift back 
    else
        propagated_wave = zeros([size(k_z), numel(distance)]);
        for n = 1:numel(distance)
            propagated_wave(:,:,n) = fftshift(ifft2((ft_wave_2d .* exp(1i * mask .* k_z .* distance(n)) .* mask))); % shift back 
        end
    end
    
    if plot_angular_spectrum
        plot_wavefunction_2d(ft_wave_2d, k_x, k_y)
    end
end

function [u, v, focal_plane_wave_2d, space_width] = get_focal_plane_wavefunction(sampled_wave_2d, adjustment, space_width, number_of_samples, k, plot_focal_plane)
    focal_plane_wave_2d = fftshift(fft2(ifftshift(sampled_wave_2d))); % transform to focal plane
    integer_list = floor(-number_of_samples/2+0.5):floor(number_of_samples/2-0.5);
    [u, v] = meshgrid(integer_list * 2 * pi / space_width / k * 8e-3); % new positions in the focal plane   
    space_width = (max(u(:)) - min(u(:))) * adjustment;
    if plot_focal_plane
        plot_wavefunction_2d(focal_plane_wave_2d, u, v)
    end
end

function xy_xz_plot_psf(propagated_wave_2d, x, y, z_list,z_plane_no)
    figure(); 
    subplot(1,2,1)
    h = pcolor(x, y, abs(propagated_wave_2d(:,:,z_plane_no)));
    set(h,'EdgeColor','none')
    axis equal

    subplot(1,2,2)
    [zz, xx] = meshgrid(z_list, x(500,:));
    h = pcolor(xx, zz, abs(squeeze(propagated_wave_2d(:,round(size(x,2)/2),:))));
    set(h,'EdgeColor','none')
    axis equal
end
