function [ M ] = lens( )

syms l1 l2 r1 r2 f x k

% l1 to f1 to f2 to l2

P1 = [1 l1; 0 1];
Q1 = [1 0; -1/r1 1];
PL = [1 0; 0 1];
Q2 = [1 0; 1/r2 1];
P2 = [1 l2; 0 1];

M2 = [1 0; -1/f 1];
M = P1 * M2 * P2;

M = M * [x;k]
end

