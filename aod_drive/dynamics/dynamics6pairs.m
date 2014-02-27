function [ sols ] = dynamics6pairs()

syms a1 a2 a3 a4 a5 a6

a = transpose([a1 a2 a3 a4 a5 a6]);
C = sym(cos(2*pi/3));
S = sym(sin(2*pi/3));

S2C = C*S*S;
C2S = S*C*C;
C3 = C.^3;
S3 = S.^3;

% go [1 0], -[1 0], [C S], -[C S], [C -S], -[C -S]
W1 = [0, 0, C2S, -C2S, -C2S, C2S];
W2 = [0, 0, S2C, -S2C, S2C, -S2C];
W3 = [0, 0, S3, -S3, -S3, S3];
W4 = [1, -1, C3, -C3, C3, -C3];

W = [W1;W2;W3;W4];
eqs = W*a;

sols = solve(eqs,'a1','a2','a3','a4','a5','a6');
sols = [sols.a1 sols.a2 sols.a3 sols.a4 sols.a5 sols.a6];

end

