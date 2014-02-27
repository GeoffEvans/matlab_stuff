function [ Y ] = SwapMatrixHalves( X )
% function [ Y ] = SwapMatrixHalves( X )
% Matrix passed must be of even dimension
midVal1 = size(X,1) / 2;
midVal2 = size(X,2) / 2;
Y = [X(midVal1+1:end, 1:end); X(1:midVal1, 1:end)];
Y = [Y(1:end, midVal2+1:end) Y(1:end, 1:midVal2)];
end

