function [ stretchedMat ] = stretch( inputMat, multiple )
% returns matrix with the elements duplicated 'multiple' times in the 2nd index

originalHeight = size(inputMat,1);
originalWidth = size(inputMat,2);
originalDepth = size(inputMat,3);
stretchedMat = repmat(inputMat,[multiple,1,1]); % insert new elements
stretchedMat = reshape(stretchedMat,[originalHeight,originalWidth*multiple,originalDepth]);

end

