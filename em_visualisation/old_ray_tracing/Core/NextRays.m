function [ nextRays ] = NextRays( rays, refIndexIn, interface )

% New rays start where previous stop at interface.
starts = rays.stops;

% Remember, interface shape is defined as x = f(y) NOT y = f(x).
% Grad of normal is neg recip of grad of tangent. 
shape = interface.shape;
d = 0.000000000001;
gradientsOfNormals = (shape(starts(2,:)) - shape(starts(2,:) + d)) / d;
anglesOfNormals = atan(gradientsOfNormals);

% find the angle between the ray and the normal
% note: if this is + then	grad(normal) < grad(raysIn/Out)
anglesIn = atan(rays.gradients) - anglesOfNormals;

% % Plot normals for debugging
line([starts(1,:) - 2; starts(1,:) + 2], [starts(2,:) - 2 .* gradientsOfNormals; starts(2,:) + 2 .* gradientsOfNormals], 'color', [0.3 0.3 0.3]);

% use Snell (for now) to calculate the angle between the normal and new ray
refIndexOut = interface.refractiveIndex;
anglesOut = asin( (refIndexIn/refIndexOut) * sin(anglesIn) );

% from the angle, calculate the gradient of the new ray
gradients = tan(anglesOut + anglesOfNormals);

% careful to round tiny gradients to zero so they get picked up
gradients(abs(gradients) < 0.0001) = 0;

nextRays = Rays(starts, gradients, []);

end

