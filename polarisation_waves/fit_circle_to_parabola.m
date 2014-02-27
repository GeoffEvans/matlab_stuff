function fit_circle_to_parabola()

x = -3:0.01:3;
y = x.^2 / 2;

subplot(1,2,1);
hold on;
grid on;
axis equal;
plot(x,y,'k','LineWidth',2);

theta = 0:0.01:2*pi;
xCircle = -0.8:0.4:0.8;
a = -xCircle.^3;
b = 1 + 3/2 * xCircle.^2;
colourMap = hsv(length(xCircle));
for k = 1:length(xCircle)
    plot(a(k) + Radius(xCircle(k)).*cos(theta), b(k) + Radius(xCircle(k)).*sin(theta), 'color', colourMap(k,:));
    line([xCircle(k) a(k)], [xCircle(k).^2/2 b(k)], 'color', colourMap(k,:));
end
xlabel('fraction of minimum radius');
ylabel('fraction of minimum radius');
hold off;

subplot(1,2,2);
plot(x,Radius(x));
grid on;
xlabel('parabola offset as fraction of minimum radius');
ylabel('focal length scale factor');

    function r = Radius(x)
        r = (1 + x.^2).^1.5;
    end
end

