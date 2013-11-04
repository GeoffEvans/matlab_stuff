function Test()

theta = 0:0.01:50*pi;

xArrayFun = Mapping(theta);
xForLoopArray = ForloopArray(theta);
xForLoopCell = ForloopCell(theta);
xPreMultiply = PreMultiplied(theta);
xForAll = ForAll(theta);
xTest = test(theta);
end

function f = fun(theta)
M = [cos(theta) sin(theta); -sin(theta) cos(theta)];
N = [theta.^2 -1; 5*theta 0];
f = M * N;
end

function f = Mapping(theta)
f = arrayfun(@fun, theta,'uniformoutput', false);
end

function f = ForloopArray(theta)
L = length(theta);
f = zeros(2,2,L);
for k = 1:L
    f(:,:,k) = fun(theta(k));
end
end

function f = ForloopCell(theta)
L = length(theta);
f = cell(1,L);
for k = 1:L
    f{k} = fun(theta(k));
end
end

function f = PreMultiplied(theta)
L = length(theta);
f = zeros(2,2,L);
f11 = cos(theta) .* theta.^2 + sin(theta) .* theta * 5;
f12 = - cos(theta);
f21 = - sin(theta) .* theta.^2 + cos(theta) .* theta * 5;
f22 = sin(theta);
f(1,1,:) = f11;
f(1,2,:) = f12;
f(2,1,:) = f21;
f(2,2,:) = f22;
end

function f = ForAll(theta)
L = length(theta);
f = zeros(2,2,L);
for k = 1:L
    f11 = cos(theta(k)) .* theta(k).^2 + sin(theta(k)) .* theta(k) * 5;
    f12 = - cos(theta(k));
    f21 = - sin(theta(k)) .* theta(k).^2 + cos(theta(k)) .* theta(k) * 5;
    f22 = sin(theta(k));
    f(1,1,k) = f11;
    f(1,2,k) = f12;
    f(2,1,k) = f21;
    f(2,2,k) = f22;
end
end