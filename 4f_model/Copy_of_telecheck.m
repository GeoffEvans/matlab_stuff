syms d1 d2 d3 d4 d5 f1 f2 f3 f4

f1 = 150;
f2 = 200;

% d1 = 96;
% d2 = 96;
% d3 = 124;
% d4 = 162;
% d5 = 153.67;

l1 = [1 0; -1/f1 1];
l2 = [1 0; -1/f2 1];

k1 = [1 d1; 0 1];
k2 = [1 d2; 0 1];
k3 = [1 d3; 0 1];

m = k1 * l1 * k2 * l2 * k3;
M = m + [f1/f2, 0; 0 f2/f1];
L = 682;

s = solve([M(:); L-d1-d2-d3], 'd1','d2','d3');
d = double([s.d1,s.d2,s.d3])
x=1