function sols = matrixdrivepairs()

syms c1 c2 c3 c4 l1 l2 l3 l4 C S vx vy

I2 = eye(2);
I4 = eye(4);
R13 = [1 0; 0 0];
R24 = [0 0; 0 1];


P1 = [I2 I2*l1;zeros(2) I2];
P2 = [I2 I2*l2;zeros(2) I2];
P3 = [I2 I2*l3;zeros(2) I2];
P4 = [I2 I2*l4;zeros(2) I2];

Q1 = [I2 zeros(2); -c1*R13 I2];
Q2 = [I2 zeros(2); -c2*R24 I2];
Q3 = [I2 zeros(2); -c3*R13 I2];
Q4 = [I2 zeros(2); -c4*R24 I2];

M = P4*Q4*P3*Q3*P2*Q2*P1*Q1;
%  M(1:2,1:2) == 0

%go 1 -1 1 -1
D13 = P4*Q4*P3*(Q3*P2*Q2*P1*c1-I4*c3);
% D1(1,3) == D1(2,3) == 0
D24 = P4*(Q4*P3*Q3*P2*c2-I4*c4);
% C D1(1,3) + S D1(1,4) == C D1(2,3) + S D1(2,4) == 0

D24q = D24(1:2,3:4) * [0;1];
D13q = D13(1:2,3:4) * [1;0];

D = D24q + D13q - [vx; vy];

eqs = [M(1:2,1:2) D];
eqs = eqs(:);
s = solve(eqs,'c1','c2','c3','c4');
sols = [s.c1; s.c2;s.c3;s.c4];

end

