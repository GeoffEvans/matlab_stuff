function [ Y ] = SwapArrayHalves( X )
% function [ Y ] = SwapArrayHalves( X )
% Array passed must be of even dimension
midVal = size(X,2) / 2;
Y = [X(midVal+1:end) X(1:midVal)];

end

