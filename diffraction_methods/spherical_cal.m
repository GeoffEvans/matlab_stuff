n = 9;
number_of_samples = 2^n - 1;
[u, v] = meshgrid(5e-5 .* linspace(-0.5, 0.5, number_of_samples));
dist = -150e-6;
k_0 = 2*pi/800e-9;

sigma = 0.1e-6;
wavefunc = @(u,v) exp(1i * k_0 * sqrt(dist.^2 + u.^2 + v.^2) - (u.^2 + v.^2)/(2*sigma.^2));
sampled_wave_2d = wavefunc(u, v);
figure; pcolor(angle(sampled_wave_2d));
focal_plane_wave_2d = fftshift(fft2(ifftshift(sampled_wave_2d)));
figure; pcolor(angle(focal_plane_wave_2d));
figure; pcolor(abs(focal_plane_wave_2d));
f_tube = 160e-3;
obj_mag = 20;
k = u * k_0 / f_tube * obj_mag;

scaling = 0.8;
integer_list = floor(-number_of_samples/2+0.5):floor(number_of_samples/2-0.5);
space_width = (max(k(:)) - min(k(:)));
x_discrete = integer_list * 2 * pi / space_width / scaling;

[x, y] = meshgrid(x_discrete);
f = unwrap(angle(focal_plane_wave_2d(2^(n-1),:)));

figure; plot(x(1,:), f, 'k');
%p = polyfit(x(1,87:169), f(87:169), 2);
%hold on 
%plot(x(1,87:169), polyval(p, x(1,87:169)), 'r--');