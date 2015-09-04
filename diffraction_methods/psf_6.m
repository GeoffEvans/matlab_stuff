function psf_6()
    time = 1e-6;
    w2 = 0;
    w4 = -2;
    ws = 4;
    a1 = [0, w2, 0, w4, 0];
    a2 = [0, w2, 0, w4, 0];    
    a3 = [0, w2, 0, w4, 0];
    a4 = [0, w2, 0, w4, 0];    
    a5 = [0, w2, 0, w4, 0];    
    a6 = [0, w2, 0, w4, 0];    
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {a1, a2, a3, a4, a5, a6};
    
    a = aol_fft();
    a.z_list = a.z_list + 0e-6;
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    a.xy_xz_plot_psf(propagated_wave, x, y, 0.6);
    res = a.get_psf_dimensions(propagated_wave, x, y, a.z_list)
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