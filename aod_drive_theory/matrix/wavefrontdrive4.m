function sols = matrixdrive6cyclic()

syms c1 c2 c3 c4 c5 c6 f vx vy L t

C = 0;
S = 1;
D1 = [1;0];
D2 = [C;S];
R13 = [1 0; 0 0];
R24 = [C*C S*C; C*S S*S];

k = [0;0];
x1 = [0;0];
[x2,k1] = prop_ray(x1,k,c1,L,D1,t);
[x3,k2] = prop_ray(x2,k1,c2,L,D2,t);
[x4,k3] = prop_ray(x3,k2,c3,L,-D1,t);
[xF,~] = prop_ray(x4,k3,c4,f,-D2,t);
x = diff(xF, 't'); % checked correct

w1 = propagate(0  + c1*R13, L);
w2 = propagate(w1 + c2*R24, L);
w3 = propagate(w2 + c3*R13, L);
w4 = propagate(w3 + c4*R24, 0);
m = simplify(w4) - eye(2)/f; % checked correct

eqs = [m x];
eqs = simplify(eqs(:));
sols = solve(eqs,'c1','c2','c3','c4');
end

function wp = propagate(w, d)
    A = eye(2) - d*w;
    wp = (w - d*det(w)*eye(2)) / det(A);
end

function [x2, k2] = prop_ray(x1, k1, w, d, D, t) % from entrance to entrance
    k2 = k1 - D * w * (dot(D,x1) - t);
    x2 = x1 + k2 * d;
end