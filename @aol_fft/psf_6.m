function psf_6(k, number_of_samples, adjustment1, adjustment2, z_list, z_plane_no, do_plot)
    waves = cell(6,1);
    scaling = 20e6;
    waves{1} = @(t) t.^2 * scaling;
    waves{2} = @(t) t.^2 * scaling;    
    waves{3} = @(t) t.^2 * scaling;
    waves{4} = @(t) t.^2 * scaling;    
    waves{5} = @(t) t.^2 * scaling;    
    waves{6} = @(t) t.^2 * scaling;    
    time = 0e-6;

    a = aol_fft();
    [sampled_wave_2d, space_width] = get_sampled_wavefunction6(waves, time, a.adjustment1, a.number_of_samples, a.k, a.do_plot.input, a.do_plot.angular_spec);  
    a.calculate_and_plot(sampled_wave_2d, space_width);
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