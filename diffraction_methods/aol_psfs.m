function aol_psfs()
    adjustment1 = 5e2; % ratio of k to max(kx)
    adjustment2 = 1e0;
    number_of_samples = 2.^10 - 1; % computational speed vs accuracy
    k = 2*pi/920e-9;
    z_list = linspace(0,100,99)*1e-6;
    z_plane_no = round(size(z_list,2)/2); % which z plane is plotted for the xy psf 
    do_plot = struct('input', 0, 'focal_plane', 0, 'angular_spec', 0);
    
    aol4_psfs(k, number_of_samples, adjustment1, adjustment2, z_list, z_plane_no, do_plot);
    %aol6_psfs(k, number_of_samples, adjustment1, adjustment2, z_list, z_plane_no, do_plot);
end

function aol4_psfs(k, number_of_samples, adjustment1, adjustment2, z_list, z_plane_no, do_plot)
    waves = struct();
    waves.x1 = [0, 10, 0, 0, 0];
    waves.y1 = [0, 10, 0, 0, 0];    
    waves.x2 = [0, 10, 1, 0, 0];
    waves.y2 = [0, 10, 0.1, 0, 0];    
    waves.focus = 0;
    waves.spherical = 0;
    waves.ac_angle = 0; % walk off angle in degrees
    time = 0e-6;
    
    [sampled_wave_2d, space_width] = get_sampled_wavefunction4(waves, time, adjustment1, number_of_samples, k, do_plot.input, do_plot.angular_spec);  
    calculate_and_plot(sampled_wave_2d, space_width, number_of_samples, k, adjustment2, z_list, z_plane_no, do_plot);
end

function aol6_psfs(k, number_of_samples, adjustment1, adjustment2, z_list, z_plane_no, do_plot)
    waves = cell(6,1);
    scaling = 20e6;
    waves{1} = @(t) t.^2 * scaling;
    waves{2} = @(t) t.^2 * scaling;    
    waves{3} = @(t) t.^2 * scaling;
    waves{4} = @(t) t.^2 * scaling;    
    waves{5} = @(t) t.^2 * scaling;    
    waves{6} = @(t) t.^2 * scaling;    
    time = 0e-6;

    [sampled_wave_2d, space_width] = get_sampled_wavefunction6(waves, time, adjustment1, number_of_samples, k, do_plot.input, do_plot.angular_spec);  
    calculate_and_plot(sampled_wave_2d, space_width, number_of_samples, k, adjustment2, z_list, z_plane_no, do_plot);
end

function calculate_and_plot(sampled_wave_2d, space_width, number_of_samples, k, adjustment2, z_list, z_plane_no, do_plot)
    [x, y, focal_plane_wave_2d, space_width] = get_focal_plane_wavefunction(sampled_wave_2d, adjustment2, space_width, number_of_samples, k, do_plot.focal_plane);
    propagated_wave_2d = propagate_wave(focal_plane_wave_2d, z_list, space_width, k, number_of_samples, do_plot.angular_spec);
    xy_xz_plot_psf(propagated_wave_2d, x, y, z_list, z_plane_no);
end

function plot_wavefunction_2d(wave_function_2d, x, y, labels)
    figure()
    subplot(1,2,1)
    hh = pcolor(x,y,abs(wave_function_2d));
    xlabel(labels{1})
    ylabel(labels{2})
    set(hh,'EdgeColor','none')
    %axis equal
    subplot(1,2,2)
    hh = pcolor(x,y,angle(wave_function_2d));
    xlabel(labels{1})
    ylabel(labels{2})
    set(hh,'EdgeColor','none')
    %axis equal
end

function [sampled_wave, space_width] = get_sampled_wavefunction4(waves, t, adjustment1, number_of_samples, k, plot_input, plot_angular_spectrum)     
    space_width = pi * number_of_samples / k * adjustment1;
    samples = linspace(-1/2, 1/2, number_of_samples) * space_width;
    [x, y] = meshgrid(samples);
    
    % take aperture to be 12mm because of scaling by relay
    scale_factors = 2*pi ./ (6e-3).^(1:5) / k;
    phase_x1 = k * (waves.x1 .* scale_factors);
    phase_y1 = k * (waves.y1 .* scale_factors);
    phase_x2 = k * (waves.x2 .* scale_factors);
    phase_y2 = k * (waves.y2 .* scale_factors);   
 
    envelope = exp(-(x.^2 + y.^2)/2/(6e-3).^2);% .* (sqrt(x.^2 + y.^2) < 5e-3);
    V = 613;
    labels = {'x', 'y'};
    sampled_wave = envelope... 
     .* exp(1i * (phase_x1(1) .* (x-V*t) + phase_x1(2) .* (x-V*t).^2 + phase_x1(3) .* (x-V*t).^3 + phase_x1(4) .* (x-V*t).^4));
    if plot_input
        plot_wavefunction_2d(sampled_wave, x, y, labels)
    end
    sampled_wave = propagate_wave(sampled_wave, 0.05, space_width, k, number_of_samples, plot_angular_spectrum)...
     .* exp(1i * (phase_y1(1) .* (y-V*t) + phase_y1(2) .* (y-V*t).^2 + phase_y1(3) .* (y-V*t).^3 + phase_y1(4) .* (y-V*t).^4)); 
    if plot_input
        plot_wavefunction_2d(sampled_wave, x, y, labels)
    end
    sampled_wave = propagate_wave(sampled_wave, 0.05, space_width, k, number_of_samples, plot_angular_spectrum)...
     .* exp(1i * (phase_x2(1) .* (-x-V*t) + phase_x2(2) .* (-x-V*t).^2 + phase_x2(3) .* (-x-V*t).^3 + phase_x2(4) .* (-x-V*t).^4));
    if plot_input
        plot_wavefunction_2d(sampled_wave, x, y, labels)
    end
    sampled_wave = propagate_wave(sampled_wave, 0.05, space_width, k, number_of_samples, plot_angular_spectrum)...
     .* exp(1i * (phase_y2(1) .* (-y-V*t) + phase_y2(2) .* (-y-V*t).^2 + phase_y2(3) .* (-y-V*t).^3 + phase_y2(4) .* (-y-V*t).^4)); 
    if plot_input
        plot_wavefunction_2d(sampled_wave, x, y, labels)
    end
    sampled_wave = sampled_wave .* exp(1i * (...  
        (x.^2 + y.^2) * waves.focus * 2*pi ./ (6e-3).^2 ./ 2 ...  
        + tand(waves.ac_angle) * 0.5 .* (x.^3 + y.^3) .* (2 * (phase_x1(2) + phase_x2(2)) / k).^2  * k/(2*pi)...  
        + (x.^2 + y.^2).^2 * waves.spherical * 2*pi ./ (6e-3).^4 ./ 4)); 
    if plot_input
        plot_wavefunction_2d(sampled_wave, x, y, labels)
    end
end

function [sampled_wave, space_width] = get_sampled_wavefunction6(waves, t, adjustment1, number_of_samples, k, plot_input, plot_angular_spectrum)     
    space_width = pi * number_of_samples / k * adjustment1;
    samples = linspace(-1/2, 1/2, number_of_samples) * space_width;
    [x, y] = meshgrid(samples);
    
    V = 613;
    labels = {'x', 'y'};
    angles = linspace(0,5*pi/3,6);
    wave_func = waves{1};
    sampled_wave = exp(-(x.^2 + y.^2)/2/(6e-3).^2) ...
        .* exp(1i * 2 * pi * wave_func(t - (cos(angles(1)) .* x + sin(angles(1)) .* y))/V); 
    if plot_input
        plot_wavefunction_2d(sampled_wave, x, y, labels)
    end
    for n = 2:6
        wave_func = waves{n};
        sampled_wave = propagate_wave(sampled_wave, 0.05, space_width, k, number_of_samples, plot_angular_spectrum)...
         .* exp(1i * 2 * pi * wave_func(t - (cos(angles(n)) .* x + sin(angles(n)) .* y))/V); 
        if plot_input
            plot_wavefunction_2d(sampled_wave, x, y, labels)
        end
    end
end

function propagated_wave = propagate_wave(wave_2d, distance, space_width, k, number_of_samples, plot_angular_spectrum)
    count = [0:floor((number_of_samples - 1)/2) floor(-(number_of_samples - 1)/2):-1];
    [k_x, k_y] = meshgrid(count * 2 * pi / space_width);
    ft_wave_2d = fft2(ifftshift(wave_2d)); % recentre wave with fftshift
       
    %display(max(k_x(:))/k);
    %display((k_y(2)-k_y(1))/k)
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
        plot_wavefunction_2d(fftshift(ft_wave_2d), fftshift(k_x), fftshift(k_y), {'kx', 'ky'})
    end
end

function [u, v, focal_plane_wave_2d, space_width] = get_focal_plane_wavefunction(sampled_wave_2d, adjustment2, space_width, number_of_samples, k, plot_focal_plane)
    focal_plane_wave_2d = fftshift(fft2(ifftshift(sampled_wave_2d))); % transform to focal plane
    integer_list = floor(-number_of_samples/2+0.5):floor(number_of_samples/2-0.5);
    [u, v] = meshgrid(integer_list * 2 * pi / space_width / k * 8e-3); % new positions in the focal plane   
    space_width = (max(u(:)) - min(u(:))) * adjustment2;
    if plot_focal_plane
        plot_wavefunction_2d(focal_plane_wave_2d, u, v, {'x_focus', 'y_focus'})
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
%    axis equal
end
