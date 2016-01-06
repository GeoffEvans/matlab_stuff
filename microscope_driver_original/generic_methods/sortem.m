function [sortedVectors,descendingArrayOfDiagonalElements]=sortem(P,D)
% this function takes in two matrices P and D, presumably the output 
% from Matlab's eig function, and then sorts the columns of P to 
% match the sorted columns of D (going from largest to smallest)
% 
% EXAMPLE: 
% 
% D =
%    -90     0     0
%      0   -30     0
%      0     0   -60
% P =
%      1     2     3
%      1     2     3
%      1     2     3
% 
% [P,D]=sortem(P,D)
% P =
%      2     3     1
%      2     3     1
%      2     3     1
% D =
%    -30     0     0
%      0   -60     0
%      0     0   -90

[descendingArrayOfDiagonalElements, indicesOfSortedArray]=sort(diag(D),'descend');
sortedVectors=P(:,indicesOfSortedArray); % arrange the columns in this order

