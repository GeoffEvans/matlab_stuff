function [ intersections, exists ] = Intersections( rays, interface )

m = rays.gradients';
x0 = rays.starts(1,:)';
y0 = rays.starts(2,:)';
c = y0 - m .* x0; 
f = interface.shape;
p = interface.position;
zeroFunctions = @(Y) m * (f(Y) + p) + c * ones(1, length(Y)) - ones(size(m,1), 1) * Y;
% zeroFunctions(Y) produces a matrix of values
% / f1(Y1)  f1(Y2)  ...  f1(Yn) \
% | f2(Y1)  f2(Y2)  ...  f2(Yn) |
% . ...     ...     ...  ...    .  
% \ fn(Y1)  ...     ...  fn(Yn) /

global maxY minY stepY;
Y = minY:stepY:maxY;
searchGrid = abs(zeroFunctions(Y)); 

% Small values => near intersection: search rows for minimums
[minVal, minIndexes] = min(searchGrid, [], 2);

yIntersections = (minIndexes - 1) .* stepY + minY;
xIntersections = (yIntersections - y0) ./ m + x0;
xIntersections(m == 0) = f(y0(m==0)) + p;
intersections = [xIntersections'; yIntersections'];

% If an intersection exists, the maximum value of the zeroFunction is
% stepY * (1 - gradientRay/gradientInterface) 
% so for 'relevant' rays and an interface which is approximately straight 
% between steps the minVal will be <= stepY.
exists = (minVal <= stepY)';

end

