function [propagated_wave_2d, samples] = propagate_wave_2d( wave_function_2d, distance, plot1, plot2, plot3 )

k = 2*pi/920e-9;
number_of_samples = 1000;
space_width = pi * number_of_samples / k;
display(space_width/number_of_samples)
samples = linspace(-1/2, 1/2, number_of_samples) * space_width;
[samples_x, samples_y] = meshgrid(samples);
sampled_wave_2d = wave_function_2d(samples_x, samples_y);
if plot1
    figure(); 
    h = pcolor(samples_x, samples_y, abs(sampled_wave_2d));
    set(h,'EdgeColor','none')
    axis equal
end

count = [0:floor((number_of_samples - 1)/2) floor(-(number_of_samples - 1)/2):-1];
[k_x, k_y] = meshgrid(count * 2 * pi / space_width);
ft_wave_2d = fft2(fftshift(sampled_wave_2d)); % recentre wave with fftshift
if plot2
    figure(); 
    h = pcolor(fftshift(k_x), fftshift(k_y), abs(fftshift(ft_wave_2d)));
    set(h,'EdgeColor','none')
    axis equal
end

k_z = sqrt(k.^2 - k_x.^2 - k_y.^2);
propagated_wave_2d = fftshift(ifft2((ft_wave_2d .* exp(1i * k_z * distance))));
if plot3
    figure(); 
    h = pcolor(samples_x, samples_y, abs(propagated_wave_2d));
    set(h,'EdgeColor','none')
    axis equal
end

end

