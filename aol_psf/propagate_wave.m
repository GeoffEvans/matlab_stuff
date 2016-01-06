function [propagated_wave, samples] = propagate_wave( wave_function, distance )

k = 10;
number_of_samples = 1000;
space_width = pi * number_of_samples / k;
display(space_width/number_of_samples)
samples = linspace(-1/2, 1/2, number_of_samples) * space_width;
sampled_wave = wave_function(samples);
figure(); plot(samples, sampled_wave);

count = [0:floor((number_of_samples - 1)/2) floor(-(number_of_samples - 1)/2):-1];
k_x = count * 2 * pi / space_width;
ft_wave = fft(fftshift(sampled_wave)); % recentre wave with fftshift
figure(); plot(fftshift(k_x), fftshift(ft_wave));

k_z = sqrt(k.^2 - k_x.^2);
propagated_wave = fftshift(ifft((ft_wave .* exp(1i * k_z * distance))));
figure(); plot(samples, abs(propagated_wave));

end

