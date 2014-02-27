function [Y] = Triangle( X )
% function [Y] = Triangle( X )
% returns a triangle between -1/2 and 1/2 with max height 1.

Window = abs(X) <= 1/2;
Y = (1 - 2*abs(X)) .* Window; 
end

