function ellipse_grid(x, y, x_fwhm, y_fwhm)

t = linspace(0, 2*pi);
figure; hold on;
for n = 1:numel(x)
    plot(x(n) + x_fwhm(n)/2 * cos(t), y(n) + y_fwhm(n)/2 * sin(t));
end
end

