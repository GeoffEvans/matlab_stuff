function [Y] = Boxcar( X )
%function [Y] = Boxcar( X )

Y = abs(X) <= 1/2;
end

