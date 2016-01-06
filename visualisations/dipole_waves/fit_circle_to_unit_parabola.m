function fit_circle_to_parabola()

x = -2:0.01:2;
y = x.^2;

subplot(1,2,1);
hold on;
grid on;
axis equal;
plot(x,y,'k','LineWidth',2);
r = 0.5;
theta = 0:0.01:2*pi;
plot(r.*cos(theta),0.5+r.*sin(theta));
line([0 0], [0 0.5]);

xCircle = -0.8:0.4:0.8;
a = -4 * xCircle.^3;
b = 0.5 + 3 * xCircle.^2;
r = 0.5 * (1 + 4*xCircle.^2).^1.5;
colourMap = hsv(length(xCircle));
for k = 1:length(xCircle)
    plot(a(k) + r(k).*cos(theta), b(k) + r(k).*sin(theta), 'color', colourMap(k,:));
    line([xCircle(k) a(k)], [xCircle(k).^2 b(k)], 'color', colourMap(k,:));
end
hold off;

subplot(1,2,2);
x = -0.8:0.01:0.8;
r = 0.5 * (1 + 4*x.^2).^1.5;
plot(x,r);
grid on;
axis([-0.8 0.8 0 8])
xlabel('parabola offset');
ylabel('focal length scale factor');
end

