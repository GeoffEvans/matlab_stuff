function ellipse_grid(x_mesh, y_mesh, x_fwhm, y_fwhm)

x_fwhm = min(x_fwhm, 3.132);

% t = linspace(0, 2*pi);
% figure; hold on;
% for n = 1:numel(x_mesh)
%     plot(x_mesh(n) + x_fwhm(n)/2*5 * cos(t), y_mesh(n) + y_fwhm(n)/2*5 * sin(t));
% end
% axis equal
% excess = 40;
% xlim([min(x_mesh(:))-excess, max(x_mesh(:))+excess])
% ylim([min(y_mesh(:))-excess, max(y_mesh(:))+excess])
% 
% figure
rng = linspace(-250,250);
hold on
plot(rng, spline(mean(x_mesh(:,[1:5,7:12])),mean(x_fwhm(:,[1:5,7:12])),rng),'r')
%plot(x_mesh, x_fwhm, 'rx')
plot(rng, spline(mean(y_mesh([1:5,7:12],:)'),mean(y_fwhm([1:5,7:12],:)'),rng),'r')
%plot(y_mesh, y_fwhm, 'bx')

wavelength = 920e-3;
wavelength_fwhm = 7e-3 / sqrt(2);
angular_fwhm = abs(rng) ./ wavelength .* wavelength_fwhm;
plot(rng, angular_fwhm, 'k--')
ylim([0, 10])

end

