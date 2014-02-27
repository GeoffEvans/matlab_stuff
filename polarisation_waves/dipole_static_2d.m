function dipole_static_2d()
permittivity = 8.85418782e-12; % m-3 kg-1 s4 A2
p = 1; % dipole moment

thetaRange = 0:0.01*pi:2*pi;
xRange = -3:0.1:3;
zRange = -3:0.1:3;
[x, z] = GenerateAllCombinations(xRange, zRange);

theta = angle(z + 1i * x);
r = sqrt(x.^2 + z.^2);
filter = r > 0.7;

Er = 2 * p .* cos(theta) ./ (4*pi*permittivity * r.^3);
Et = p .* sin(theta) ./ (4*pi*permittivity * r.^3);

Ex = (Er .* sin(theta) - Et .* cos(theta));
Ez = (Er .* cos(theta) + Et .* sin(theta));

quiver(x(filter),z(filter),Ex(filter),Ez(filter));
hold on;

for k = 1:10
    rPol = (4 * cos(thetaRange).^2 + sin(thetaRange).^2).^(1/3)*(0.5+k/5);
    plot(rPol.*sin(thetaRange), rPol.*cos(thetaRange), 'r');
end
hold off;

end

function [r, theta] = GenerateAllCombinations(rRange,thetaRange)
    [ thetaMesh, rMesh ] = meshgrid(thetaRange, rRange);
    theta = thetaMesh(:);
    r = rMesh(:);
end

