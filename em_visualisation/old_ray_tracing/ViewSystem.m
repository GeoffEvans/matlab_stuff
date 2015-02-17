function [] = ViewSystem(interfaces, rays)
% Produces a 2D plot of an optical system

global opticalLength;
currentRefIndex = 1;
endPlane = Interface(opticalLength, @Constant, 1);
interfaces = [interfaces, endPlane];

grid;
grid minor;
hold on;
PlotInterfaces(interfaces);

for j = 1: length(interfaces)
    thisInterface = interfaces(j);
    [rays, currentRefIndex] = PropagateRays(rays, currentRefIndex, thisInterface);
end
hold off;
end

function [newRays, newIndex] = PropagateRays(rays, currentRefIndex, interface)
    [rays.stops, exists] = Intersections(rays, interface);    
    raysFiltered = Rays(rays.starts(:,exists), rays.gradients(exists), rays.stops(:,exists)); % Remove rays that don't intersect with the interface
    PlotRays(raysFiltered);
    newRays = NextRays(raysFiltered, currentRefIndex, interface);
    newIndex = interface.refractiveIndex;
end