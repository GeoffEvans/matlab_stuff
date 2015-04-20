function [ tRamp, A, B, C ] = scaleTrampABC( aP,bP,cP,tScan,systemVar )
% Check the array dimensions inorder to use the same function for all
% Modes.
A = round( aP * (2^32/systemVar.controlClockFreq)); 
%A(:,[2,3]) = A(:,[3,2]);

Bmy = round(bP .* (2^32/systemVar.controlClockFreq^2));
%Bmy(:,[2,3]) = Bmy(:,[3,2]);

% scale B to allow for steep slopes withing 16 bits
B = round(Bmy./8); 
B(:,[2,3]) = B(:,[3,2]);

C = round(cP .* 8192. * (2^32/systemVar.controlClockFreq^3));
C(:,[2,3]) = C(:,[3,2]);

tRamp = round(tScan .* systemVar.controlClockFreq);
tRamp(:,[2,3]) = tRamp(:,[3,2]);
end