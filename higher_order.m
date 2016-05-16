% define a profile, propagate a distance, predict and compute profiles
Q = 1e-9;
L = 1;
num_rays = 101;

w1 = 0;
w2 = 1/3;
w3 = 0;
w4 = 0;
w5 = 0;

w = @(x) w1*x + w2*x.^2/2 + w3*x.^3/6 + w4*x.^4/24 + w5*x.^5/120;

rays_start = linspace(-0.01, 0.01, num_rays);
rays_direction = -(w(rays_start + Q) - w(rays_start - Q)) / 2 / Q;
rays_stop = rays_start + rays_direction * L;

gamma = 1 ./ (1 - L*w2);
w_stop1 = @(x) w1*x + gamma*w2*x.^2/2 + gamma.^3*w3*x.^3/6 + gamma.^4*w4*x.^4/24 + gamma^5*w5*x.^5/120;
w_stop2 = @(x) w1*x + gamma*w2*x.^2/2 + gamma.^3*w3*x.^3/6 + gamma.^4*(w4 + 3*L*(gamma*w3^2))*x.^4/24 + gamma^5*w5*x.^5/120;
w_stop3 = @(x) w1*x + gamma*w2*x.^2/2 + gamma.^3*w3*x.^3/6 + gamma.^4*(w4 + 3*L*gamma*w3^2)*x.^4/24 + gamma^5*(w5 + 5*gamma*L*w3*(2*w4+3*L*gamma*w3^2))*x.^5/120;
w_stop3 = @(x) w1*x + gamma*w2*x.^2/2 + gamma.^3*w3*x.^3/6 + gamma.^4*(w4 + 3*L*gamma*w3^2)*x.^4/24 + gamma^5*(w5 + 5*gamma*L*w3*(2*w4+3*L*gamma*w3^2-6*w2^3))*x.^5/120;

figure()
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
subplot(1,2,1)
plot([rays_start; rays_stop], [zeros(1, num_rays); L * ones(1, num_rays)]);
subplot(1,2,2)
hold on
plot(rays_start, w(rays_start), 'r-')
plot(rays_stop, w(rays_start), 'o')
plot(rays_start, w_stop3(rays_start), 'm')
plot(rays_start, w_stop2(rays_start), 'b-.')
plot(rays_start, w_stop1(rays_start), 'g--')
ylim([0, 4e-5])

hold off