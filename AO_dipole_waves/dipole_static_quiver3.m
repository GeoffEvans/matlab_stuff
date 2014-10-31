function dipole_static()
permittivity = 8.85418782e-12; % m-3 kg-1 s4 A2
p = 1; % dipole moment

rRange = 0.1:0.01:1;
thetaRange = 0:0.01*pi:pi;
phiRange = 0:0.1*pi:2*pi;
[r, theta, phi] = GenerateAllCombinations(rRange,thetaRange,phiRange);

Er = 2 * p .* cos(theta) ./ (4*pi*permittivity * r.^3);
Et = p .* sin(theta) ./ (4*pi*permittivity * r.^3);

Ex = ( Er .* sin(theta) - Et .* cos(theta) ) .* cos(phi);
Ey = ( Er .* sin(theta) - Et .* cos(theta) ) .* sin(phi);
Ez = ( Er .* cos(theta) + Et .* sin(theta) );

x = sin(theta) .* cos(phi);
y = sin(theta) .* sin(phi);
z = cos(theta);

quiver3(x,y,z,Ex,Ey,Ez);

end

function [r, theta, phi] = GenerateAllCombinations(rRange,thetaRange,phiRange)
    [ thetaMeshTemp, rMeshTemp ] = meshgrid(thetaRange, rRange);
    [ phiMesh, thetaMesh ] = meshgrid(phiRange, thetaMeshTemp(:));
    [ ~, rMesh ] = meshgrid(phiRange, rMeshTemp(:));
    theta = thetaMesh(:);
    phi = phiMesh(:);
    r = rMesh(:);
end