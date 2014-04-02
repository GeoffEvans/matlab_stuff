function [ x,y,wavefront ] = find_aberrations_on_rays( rayBundle )
% Only run this on functions with single time, single drive, single perturbation

% CURRENTLY ASSUMES THAT FOCUS IS ON AXIS
xFocal = 0;
yFocal = 0;

focalLength = rayBundle.zFocusPredicted;
rayWavevector = rayBundle.k;
rayPositionAol = rayBundle.GetXyzLeavingAol();

rayPositionFocalPlane = propagate_ray_to_plane(rayPositionAol,rayWavevector,stretch([0;0;1],rayBundle.numOfRays), stretch([0;0;focalLength],rayBundle.numOfRays));

x = rayPositionFocalPlane(1,:,:);
y = rayPositionFocalPlane(2,:,:);

dWdx = (xFocal - x) / focalLength;
dWdy = (yFocal - y) / focalLength;

W = 

wavefront = -( sqrt(focalLength.^2 + x.^2 + y.^2) + W );

end

