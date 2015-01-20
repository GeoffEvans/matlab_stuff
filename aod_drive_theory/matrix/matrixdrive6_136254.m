function sols = matrixdrive6cyclic(order)
syms c1 c2 c3 c4 c5 c6 f vx vy

a = order;

L = 5e-2;

l1 = L;
l2 = L;
l3 = L;
l4 = L;
l5 = L;
l6 = f;

I2 = eye(2);
I4 = eye(4);
C = cos(2*pi/3);
S = sin(2*pi/3);

R = cell(6,1);
R{1} = [1 0; 0 0];
R{2} = [C*C -S*C; -C*S S*S];
R{3} = [C*C S*C; C*S S*S];
R{4} = R{1};
R{5} = R{2};
R{6} = R{3};

q = cell(6,1);
q{1} = [1;0];
q{2} = [-C;S];
q{3} = [C;S];
q{4} = -q{1};
q{5} = -q{2};
q{6} = -q{3};

P1 = [I2 I2*l1;zeros(2) I2];
P2 = [I2 I2*l2;zeros(2) I2];
P3 = [I2 I2*l3;zeros(2) I2];
P4 = [I2 I2*l4;zeros(2) I2];
P5 = [I2 I2*l5;zeros(2) I2];
P6 = [I2 I2*l6;zeros(2) I2];

Q1 = [I2 zeros(2); -c1*R{a(1)} I2];
Q2 = [I2 zeros(2); -c2*R{a(2)} I2];
Q3 = [I2 zeros(2); -c3*R{a(3)} I2];
Q4 = [I2 zeros(2); -c4*R{a(4)} I2];
Q5 = [I2 zeros(2); -c5*R{a(5)} I2];
Q6 = [I2 zeros(2); -c6*R{a(6)} I2];

M = P6*Q6*P5*Q5*P4*Q4*P3*Q3*P2*Q2*P1*Q1;

D1 = P6*Q6*P5*Q5*P4*Q4*P3*Q3*P2*Q2*P1*c1;
D2 = P6*Q6*P5*Q5*P4*Q4*P3*Q3*P2*c2;
D3 = P6*Q6*P5*Q5*P4*Q4*P3*c3;
D4 = P6*Q6*P5*Q5*P4*c4;
D5 = P6*Q6*P5*c5;
D6 = P6*c6;

D1q = D1(1:2,3:4) * q{a(1)};
D2q = D2(1:2,3:4) * q{a(2)};
D3q = D3(1:2,3:4) * q{a(3)};
D4q = D4(1:2,3:4) * q{a(4)};
D5q = D5(1:2,3:4) * q{a(5)};
D6q = D6(1:2,3:4) * q{a(6)};

D = D1q + D2q + D3q + D4q + D5q + D6q - [vx;vy];

eqs = [M(1:2,1:2) D];
eqs = eqs(:);
sols = solve(eqs,'c1','c2','c3','c4','c5','c6');
%sols = [sols.c1; sols.c2; sols.c3; sols.c4; sols.c5; sols.c6];
end

