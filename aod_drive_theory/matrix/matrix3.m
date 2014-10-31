function sF = matrix3()

syms c1 c2 c3 l1 l2 l3 a2 a3 vx vy

I2 = eye(2);
C2 = cos(a2);
S2 = sin(a2);
C3 = cos(a3);
S3 = sin(a3);
R1 = [1 0; 0 0];
R2 = [C2*C2 S2*C2; C2*S2 S2*S2];
R3 = [C3*C3 S3*C3; C3*S3 S3*S3];

P1 = [I2 I2*l1;zeros(2) I2];
P2 = [I2 I2*l2;zeros(2) I2];
P3 = [I2 I2*l3;zeros(2) I2];

Q1 = [I2 zeros(2); -c1*R1 I2];
Q2 = [I2 zeros(2); -c2*R2 I2];
Q3 = [I2 zeros(2); -c3*R3 I2];

%%%%%%%%%%%%%%%
D1 = P3*Q3*P2*Q2*P1*c1;
D2 = P3*Q3*P2*c2;
D3 = P3*c3;

D1q = D1(1:2,3:4) * [1;0];
D2q = D2(1:2,3:4) * [C2;S2];
D3q = D3(1:2,3:4) * [C3;S3];

D = D3q + D2q + D1q;
M = P3*Q3*P2*Q2*P1*Q1;
eqs = M(1:2,1:2);

s = solve([eqs(:); D(:)], 'c1', 'c2', 'c3', 'a2', 'a3');

sols = [s.c1; s.c2; s.c3]

end


