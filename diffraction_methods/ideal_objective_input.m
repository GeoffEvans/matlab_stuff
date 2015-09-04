function ideal_objective_input()
k = 2*pi/850e-9;
focal_disp = 100e-6;
space_width_out = 5e-2;
space_width_in = 30e-5;
num_samples = 2^8 - 1;
[u, v] = meshgrid(linspace(-space_width_in/2, space_width_in/2, num_samples));
r = sqrt(u.^2 + v.^2);
focal_plane_wave_2d = exp(-r.^2/2/1e-10) .* exp(1i * k * sqrt(r.^2 + focal_disp.^2));

% figure()
% subplot(1,2,1)
% pcolor(u,v,abs(focal_plane_wave_2d)); 
% axis equal
% subplot(1,2,2)
% pcolor(u,v,angle(focal_plane_wave_2d)); 
% axis equal

ideal_input = fftshift(ifft2(ifftshift(focal_plane_wave_2d)));
[x, y] = meshgrid(linspace(-space_width_out/2, space_width_out/2, num_samples)); % new positions in the focal plane

% figure()
% subplot(1,2,1)
% pcolor(x,y,abs(ideal_input)); 
% axis equal
% subplot(1,2,2)
% pcolor(x,y,angle(ideal_input)); 
% axis equal

a = x(x > 0 & y == 0)
b = angle(ideal_input(x > 0 & y == 0))
c = [0; b(2:end-1) - b(1:end-2) > pi];
plot(a(1:end-1),b(1:end-1) - 2*pi*cumsum(c))
polyfit(a(1:end-1), b(1:end-1) - 2*pi*cumsum(c), 8)
end