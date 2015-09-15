function [U,V,W] = stuff()
I = eye(3);
M = [-1, 0, 0; 0, 1, 0; 0, 0, 1];

syms a b c d e f u v u1 u2 u3 v1 v2 v3

t = [sqrt(1-u^2-v^2) u v; a b c; d e f];
q = solve([inv(t)==transpose(t), det(t)==1], a, b, c, d, e, f);
T = simplify(subs(t, {a,b,c,d,e,f}, {q.a(1), q.b(1), q.c(1), q.d(1), q.e(1), q.f(1)}));
T1 = subs(T, {u, v}, {u1, 0});
T2 = subs(T, {u, v}, {u2, v2});
T3 = subs(T, {u, v}, {u3, v3});

U = simplify(transpose(T3) * M * T3 * transpose(T2) * M * T2 * transpose(T1) * M * T1); 
V = simplify(transpose(T2) * M * T2 * transpose(T1) * M * T1); 
W = simplify(transpose(T1) * M * T1); 
%s = solve(V*[1;0;0] == -[1;0;0], u1, u2, u3, v2, v3);
%t = solve(W(1) == -1, u1, u2, v2);


%subs(simplify(subs(U, {'u1','v2','u3','v3'}, [1/sqrt(2),1/sqrt(2),0,-1/sqrt(2)])),'u2',-1/sqrt(2))
end