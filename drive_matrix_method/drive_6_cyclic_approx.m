function procd = drive_6_cyclic_approx()

syms c1 c2 c3 c4 c5 c6 f vx vy
syms L positive

L = 4.554e-2

l1 = L;
l2 = L;
l3 = L;
l4 = L;
l5 = L;
l6 = f;

I2 = sym(eye(2));
I4 = sym(eye(4));
C = cos(sym(2*pi/3));
S = sin(sym(2*pi/3));
R14 = [C*C S*C; C*S S*S];
R25 = [1 0; 0 0];
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
D14 = P6*Q6*P5*Q5*P4*(Q4*P3*Q3*P2*Q2*P1*c1-I4*c4);
D25 = P6*Q6*P5*(Q5*P4*Q4*P3*Q3*P2*c2-I4*c5);
D36 = P6*(Q6*P5*Q5*P4*Q4*P3*c3-I4*c6);

D36q = D36(1:2,3:4) * [C;-S];
D25q = D25(1:2,3:4) * -[1;0];
D14q = D14(1:2,3:4) * [C;S];

D = D36q + D14q + D25q - 0.5*[vx; -vx.*sqrt(3)];

eqs = [M(1:2,1:2) D];
eqs = (eqs(:));
sols = solve(eqs,'c1','c2','c3','c4','c5','c6');

fields = fieldnames(sols);
procd = sym(zeros(6,1));
for k  = 1:6
    eqs = sols.(fields{k});
    procd(k) = subs(eqs(2), 'z', '1/3/f'); % take second solution and set to symmetric split
end  
end

