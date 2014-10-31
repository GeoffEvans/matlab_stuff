function N = matrixdrive()

syms F1 F2 F3 l1 l2 l3 L1 L2 L3

% l1 = 0;
% l2 = 0;
% l3 = 0;
% L1 = 0;
% L2 = 0;
rt3 = sqrt(sym(3));

R1 = [1 0; 0 0];
R2 = [1 -rt3; -rt3 3]/4;
R3 = [1 rt3; rt3 3]/4;

a1 = l1/(2*F1-l1);
b1 = -2/(2*F1-l1);
c1 = l1/(2*F1-2*l1);
a2 = l2/(2*F2-l2);
b2 = -2/(2*F2-l2);
c2 = l2/(2*F2-2*l2);
a3 = l3/(2*F3-l3);
b3 = -2/(2*F3-l3);
c3 = l3/(2*F3-2*l3);

D1 = [eye(2)+a1*R1 eye(2)*l1;b1*R1 eye(2)+c1*R1];
D2 = [eye(2)+a2*R2 eye(2)*l2;b2*R2 eye(2)+c2*R2];
D3 = [eye(2)+a3*R3 eye(2)*l3;b3*R3 eye(2)+c3*R3];

P1 = [eye(2) eye(2)*L1; zeros(2) eye(2)];
P2 = [eye(2) eye(2)*L2; zeros(2) eye(2)];
P3 = [eye(2) eye(2)*L3; zeros(2) eye(2)];

M = P3 * D3 * P2 * D2 * P1 * D1;
N = simplify(M(1:2,1:2));

end

