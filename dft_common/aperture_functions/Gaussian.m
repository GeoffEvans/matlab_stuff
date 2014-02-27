function [Y] = Gaussian( X )
% function [Y] = Gaussian( X )
% returns a unit Gaussian on X centred at the origin

Y = gaussmf(X, [1, 0]);
end

