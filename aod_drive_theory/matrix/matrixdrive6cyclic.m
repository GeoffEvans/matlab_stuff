function sols = matrixdrive6cyclic()

syms c1 c2 c3 c4 c5 c6 f vx vy

L = 5e-2

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
R14 = [1 0; 0 0];
R25 = [C*C S*C; C*S S*S];
R36 = [C*C -S*C; -C*S S*S];

P1 = [I2 I2*l1;zeros(2) I2];
P2 = [I2 I2*l2;zeros(2) I2];
P3 = [I2 I2*l3;zeros(2) I2];
P4 = [I2 I2*l4;zeros(2) I2];
P5 = [I2 I2*l5;zeros(2) I2];
P6 = [I2 I2*l6;zeros(2) I2];

Q1 = [I2 zeros(2); -c1*R14 I2];
Q2 = [I2 zeros(2); -c2*R25 I2];
Q3 = [I2 zeros(2); -c3*R36 I2];
Q4 = [I2 zeros(2); -c4*R14 I2];
Q5 = [I2 zeros(2); -c5*R25 I2];
Q6 = [I2 zeros(2); -c6*R36 I2];

M = P6*Q6*P5*Q5*P4*Q4*P3*Q3*P2*Q2*P1*Q1;
%  M(1:2,1:2) == 0

% go [1 0], -[C S], [C -S], -[1 0], [C S], -[C -S]
D14 = P6*Q6*P5*Q5*P4*(Q4*P3*Q3*P2*Q2*P1*c1-I4*c4);
% D1(1,3) == D1(2,3) == 0
D25 = P6*Q6*P5*(Q5*P4*Q4*P3*Q3*P2*c2-I4*c5);
% C D1(1,3) + S D1(1,4) == C D1(2,3) + S D1(2,4) == 0
D36 = P6*(Q6*P5*Q5*P4*Q4*P3*c3-I4*c6);
% C D1(1,3) - S D1(1,4) == C D1(2,3) - S D1(2,4) == 0

D36q = D36(1:2,3:4) * [C;-S];
D25q = D25(1:2,3:4) * -[C;S];
D14q = D14(1:2,3:4) * [1;0];

D = D36q + D14q + D25q - vx*[1;1];

eqs = [M(1:2,1:2) D];
eqs = eqs(:);
sols = solve(eqs,'c1','c2','c3','c4','c5','c6');

end

