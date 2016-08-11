function drive_eqs = calc_drive_eqs()

syms c1 c2 c3 c4 c5 c6 f v % six linear chirps c1-6, the focal length f, and the focal speed v*(acoustic velcity)
s = 5e-2; % spacing between AODs
velocity_unit_vector = [1;0]; % focus moves in positive x-direction with speed v, also works for 0.5*[-1, sqrt(3)] and 0.5*[-1, -sqrt(3)]

I2 = sym(eye(2)); % set up some matrices for later use
I4 = sym(eye(4));
C = cos(sym(2*pi/3));
S = sin(sym(2*pi/3));
R14 = [C*C S*C; C*S S*S];
R25 = [1 0; 0 0];
R36 = [C*C -S*C; -C*S S*S];

Ps = [I2 I2*s;zeros(2) I2]; % propagate between AODs
Pf = [I2 I2*f;zeros(2) I2]; % propagate from last AOD to focus

Q1 = [I2 zeros(2); -c1*R14 I2]; % diffract at each AOD
Q2 = [I2 zeros(2); -c2*R25 I2];
Q3 = [I2 zeros(2); -c3*R36 I2];
Q4 = [I2 zeros(2); -c4*R14 I2];
Q5 = [I2 zeros(2); -c5*R25 I2];
Q6 = [I2 zeros(2); -c6*R36 I2];

% take acoustic directions to be in order: [1 0], -[C S], [C -S], -[1 0], [C S], -[C -S]

D14 = Pf*Q6*Ps*Q5*Ps*(Q4*Ps*Q3*Ps*Q2*Ps*c1-I4*c4);
D25 = Pf*Q6*Ps*(Q5*Ps*Q4*Ps*Q3*Ps*c2-I4*c5);
D36 = Pf*(Q6*Ps*Q5*Ps*Q4*Ps*c3-I4*c6);

D36q = D36(1:2,3:4) * [C;-S];
D25q = D25(1:2,3:4) * -[1;0];
D14q = D14(1:2,3:4) * [C;S];

D = D36q + D14q + D25q; % - v * velocity_unit_vector; % matrix to constrain for focus velocity

M = Pf*Q6*Ps*Q5*Ps*Q4*Ps*Q3*Ps*Q2*Ps*Q1; % matrix to contrain for stigmatic focus: M(1:2,1:2) == 0

eqs = [M(1:2,1:2) D];
drive_eqs = solve(eqs(:),'c1','c2','c3','c4','c5','c6');

end

