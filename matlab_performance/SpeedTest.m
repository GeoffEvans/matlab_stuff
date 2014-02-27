function [xArrayFunT,xForLoopArrayT,xForLoopCellT,xPreMultiplyT,xMexT,xParforT] = SpeedTest(count)
theta = linspace(-pi,pi,count);

t = cputime; 
xArrayFun = ArrayFun(theta);
xArrayFunT = cputime-t;
t = cputime; 
xForLoopArray = ForloopArray(theta);
xForLoopArrayT = cputime-t;
t = cputime; 
xForLoopCell = ForloopCell(theta);
xForLoopCellT = cputime-t;
t = cputime; 
for N = 1:100
xPreMultiply = PreMultiplied(theta);
end
xPreMultiplyT = cputime-t;
t = cputime; 
for N = 1:100
xMex = test(theta); % Mex
end
xMexT = cputime-t;
t = cputime; 
xParfor = DoParfor(theta);
xParforT = cputime-t;
end

function f = fun(theta)
M = [cos(theta) sin(theta); -sin(theta) cos(theta)];
N = [theta.^2 -1; 5*theta 0];
f = M * N;
end

function f = ArrayFun(theta)
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

function f = DoParfor(theta)
L = length(theta);
f = zeros(L,2,2);
parfor k = 1:L
    f(k,:,:) = fun(theta(k));
end
end

