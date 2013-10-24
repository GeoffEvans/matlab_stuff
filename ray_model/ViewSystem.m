function [] = ViewSystem(interfaces, rays)

% Produces a 2D plot of an optical system

global opticalLength maxY minY stepY;

endPlane = Interface(opticalLength, @Constant, 1);
interfaces = [interfaces, endPlane];

% Configure plot
grid;
grid minor;
hold on;
PlotInterfaces(interfaces);

% Declare initial conditions
currentRefIndex = 1;

for j = 1: length(interfaces)
    thisInterface = interfaces(j);
    [rays.stops, exists] = Intersections(rays, thisInterface);
    % Remove rays that don't intersect with the interface
    raysFiltered = Rays(rays.starts(:,exists), rays.gradients(exists), rays.stops(:,exists));
    PlotRays(raysFiltered);
    rays = NextRays(raysFiltered, currentRefIndex, thisInterface);
    currentRefIndex = thisInterface.refractiveIndex;
end

hold off;

end