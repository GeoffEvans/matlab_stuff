function [ inputRays ] = PointSourceRays( focusPoint, spread, density )
% given focus point, return a rays object

% treat as a point source - uniform angular distribution
angleBetweenRaysRadians = (pi / 180) * density;
maxAngle = (pi / 180) * spread;
angles = -maxAngle:angleBetweenRaysRadians:maxAngle;

gradients = tan(angles);

starts = [zeros(1, length(gradients)); focusPoint(2) - focusPoint(1) .* gradients];
inputRays = Rays(starts, gradients, []);

end

