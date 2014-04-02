function [xyzPlane] = propagate_ray_to_plane(xyzIn,k,unitNormalToPlane,pointOnPlane)
% Propagate from a known position and direction to the surface of a plane

% the vector from a point to anywhere on a plane dotted with the normal equals the distance to the plane from the point
vectorToPointOnPlane = pointOnPlane - xyzIn;
multiplesOfWavevector = dot(vectorToPointOnPlane,unitNormalToPlane) ./ dot(k,unitNormalToPlane);
xyzPlane = xyzIn + repmat(multiplesOfWavevector,3,1) .* k; 

end