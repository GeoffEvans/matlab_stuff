function sols = matrixdrive6pairs()

syms c1 c2 c3 c4 c5 c6 l1 l2 l3 l4 l5 l6 C S L f vx vy

l1 = L;
l2 = L;
l3 = L;
l4 = L;
l5 = L;
l6 = f;

I2 = eye(2);
I4 = eye(4);
C = sym(cos(2*pi/3));
S = sym(sin(2*pi/3));
R12 = [1 0; 0 0];
R34 = [C*C S*C; C*S S*S];
R56 = [C*C -S*C; -C*S S*S];

P1 = [I2 I2*l1;zeros(2) I2];
P2 = [I2 I2*l2;zeros(2) I2];
P3 = [I2 I2*l3;zeros(2) I2];
P4 = [I2 I2*l4;zeros(2) I2];
P5 = [I2 I2*l5;zeros(2) I2];
P6 = [I2 I2*l6;zeros(2) I2];

Q1 = [I2 zeros(2); -c1*R12 I2];
Q2 = [I2 zeros(2); -c2*R12 I2];
Q3 = [I2 zeros(2); -c3*R34 I2];
Q4 = [I2 zeros(2); -c4*R34 I2];
Q5 = [I2 zeros(2); -c5*R56 I2];
Q6 = [I2 zeros(2); -c6*R56 I2];

M = P6*Q6*P5*Q5*P4*Q4*P3*Q3*P2*Q2*P1*Q1;
%  M(1:2,1:2) == 0

% go [1 0], -[1 0], [C S], -[C S], [C -S], -[C -S]
D12 = P6*Q6*P5*Q5*P4*Q4*P3*Q3*P2*(Q2*P1*c1-I4*c2);
% D1(1,3) == D1(2,3) == 0
D34 = P6*Q6*P5*Q5*P4*(Q4*P3*c3-I4*c4);
% C D1(1,3) + S D1(1,4) == C D1(2,3) + S D1(2,4) == 0
D56 = P6*(Q6*P5*c5-I4*c6);
% C D1(1,3) - S D1(1,4) == C D1(2,3) - S D1(2,4) == 0

D56q = D56(1:2,3:4) * [C;-S];
D34q = D34(1:2,3:4) * [C;S];
D12q = D12(1:2,3:4) * [1;0];

D = D56q + D34q + D12q - [vx;vy];

eqs = [M(1:2,1:2) D];
sols = solve(eqs(:),'c1','c2','c3','c4','c5','c6');
sols = [sols.c1; sols.c2; sols.c3; sols.c4; sols.c5; sols.c6];
end

