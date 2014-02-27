function [A] = Boxcar2d(X, Y)
% function [A] = Boxcar(X, Y)

A = (abs(X) <= 1/2) & (abs(Y) <= 1/2);
end

