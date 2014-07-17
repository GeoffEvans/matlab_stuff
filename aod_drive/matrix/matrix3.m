function D = matrix3()

syms c1 c2 c3 l1 l2 l3 a1 a2 vx vy

I2 = eye(2);
C1 = cos(a1);
S1 = sin(a1);
C2 = cos(a2);
S2 = sin(a2);
R1 = [1 0; 0 0];
R2 = [C1*C1 S1*C1; C1*S1 S1*S1];
R3 = [C2*C2 S2*C2; C2*S2 S2*S2];

P1 = [I2 I2*l1;zeros(2) I2];
P2 = [I2 I2*l2;zeros(2) I2];
P3 = [I2 I2*l3;zeros(2) I2];

Q1 = [I2 zeros(2); -c1*R1 I2];
Q2 = [I2 zeros(2); -c2*R2 I2];
Q3 = [I2 zeros(2); -c3*R3 I2];

M = P3*Q3*P2*Q2*P1*Q1;
%  M(1:2,1:2) == 0

eqs = [M(1:2,1:2)];
s = solve(eqs(:),'c1','c2','c3');
sols = [s.c1; s.c2;s.c3];

c1 = sols(1);
c2 = sols(2);
c3 = sols(3);

D1 = P3*Q3*P2*Q2*P1*c1;
D2 = P3*Q3*P2*c2;
D3 = P3*c3;

D3q = D3(1:2,3:4) * [C2;S2];
D2q = D2(1:2,3:4) * [C1;S1];
D1q = D1(1:2,3:4) * [1;0];

D = D3q + D2q + D1q

end


