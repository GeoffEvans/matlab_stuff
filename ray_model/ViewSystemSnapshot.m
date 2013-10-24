% Produces a 2D plot of an optical system

global opticalLength maxY minY stepY;
opticalLength = 150;
maxY = 10;
minY = -10;
stepY = 0.01;

% Declare interfaces
curvy = Interface(50, @(x) 0.2 * cos(x), 1.2);
opticalFlat = Interface(51, @Constant, 1);
interfaces = [curvy, opticalFlat, endPlane];

% Configure plot
f = figure();
hold;
grid;
grid minor;
PlotInterfaces(interfaces);

% Declare initial conditions
currentRefIndex = 1;
rays = PlaneWaveRays(@(x) 1);

for j = 1: length(interfaces)
    thisInterface = interfaces(j);
    [rays.stops, exists] = Intersections(rays, thisInterface);
    % Remove rays that don't intersect with the interface
    raysFiltered = Rays(rays.starts(:,exists), rays.gradients(exists), rays.stops(:,exists));
    PlotRays(raysFiltered);
    rays = NextRays(raysFiltered, currentRefIndex, thisInterface);
    currentRefIndex = thisInterface.refractiveIndex;
end
