function [ extendedX, extendedY ] = set_cross_product( X, Y )

% Given two matrices, X and Y, the function creates new matrices, extendedX and
% extendedY, that contain all possible permutations of the rows of X
% and Y.

% If X is of height Lx and Y is of height Ly then the extended vectors will
% be of height Lx * Ly.

sizeY = size(Y);
sizeX = size(X);
noRowsY = sizeY(1);
noRowsX = sizeX(1);
extendedX = repmat(X, sizeY(1), 1);

noRowsExtended = sizeY(1)*noRowsX;
extendedY = zeros( noRowsExtended, sizeY(2) );
for k = 1:noRowsY
    indexArray = (1:noRowsX) + noRowsX*(k-1);
    extendedY(indexArray, :) = repmat(Y(k,:),noRowsX,1);
end

end

