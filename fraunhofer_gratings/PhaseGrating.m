
% Set up a thin phase grating with an aperture function along the x direction at
% z = 0. View intensity in the xz plane.

opticWaveLength = 1e-4;
k = 2 * pi / opticWaveLength;

surfEnd = 1e1;
surfRes = opticWaveLength * 0.1; 
surface = -surfEnd:surfRes:surfEnd;
zStart = 990;
zEnd = 1400;
zRes = 10;
xEnd = 10;
xRes = 5;
[xGrid,zGrid] = meshgrid(-xEnd:xRes:xEnd, zStart:zRes:zEnd);

ApertureFunction = @(surf) surf.^2/2/1200;
R = @(X0,X,Z) (X.^2 + Z.^2).^0.5;
CosAngle = @(X0,X,Z) Z ./ R(X0,X,Z);
GreensFunc = @(X0,X,Z) exp(1i .* k .* (X - surface)/R(X0,X,Z) .* surface) ./ R(X0,X,Z);
Integrand = @(X,Z) exp(1i .* ApertureFunction(surface)) .* GreensFunc(surface,X,Z) .* CosAngle(surface,X,Z);
SurfaceIntegral = @(X,Z) sum(Integrand(X,Z));
Huygens = @(X,Z) -1i/opticWaveLength * arrayfun(SurfaceIntegral,X,Z);
FindIntensityAtXZ = @(X,Z) log(abs(Huygens(X,Z)+1));

f = figure();
subplot(1,3,1); % Plot the aperture function;
plot(surface,ApertureFunction(surface));
subplot(1,3,[2 3]); % Plot the intensity in the plane
contourf(xGrid,zGrid,FindIntensityAtXZ(xGrid,zGrid),'EdgeColor','none');
colorbar;