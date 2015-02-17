function v = matrix3()

syms c1 c2 c3 l1 l2 l3 L f vx vy x y F

l1 = 5e-2;
l2 = 5e-2;
l3 = f;

I2 = eye(2);
C = cos(pi/3);
S = sin(pi/3);
R1 = [1 0; 0 0];
R2 = [C*C S*C; C*S S*S];
R3 = [C*C -S*C; -C*S S*S];

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
D2q = D2(1:2,3:4) * [C;S];
D3q = D3(1:2,3:4) * [-C;S];

D = D3q + D2q + D1q;
M = simplify(P3*Q3*P2*Q2*P1*Q1);
eqs = M(1:2,1:2);%-M(1:2,3:4)/F;

s = solve([eqs(:)], 'c1', 'c2', 'c3');

sols = [s.c1; s.c2; s.c3]

v = subs(D, {c1,c2,c3}, {sols(1), sols(2), sols(3)});
simplify(subs(M, {c1,c2,c3}, {sols(1), sols(2), sols(3)}));
end


