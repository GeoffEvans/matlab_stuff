function doit()

options = optimoptions(@fminunc,'GradObj','on');
[x, fval] = fminunc(@min_fun, [2/3,0,2/3,0,2/3,0], options)

end

function [val, grad] = min_fun(array5)
p = stuff();
q = subs(p, {'u1','u2','u3','v1','v2','v3'}, array5);
val = norm(double(q) - [-1, 0, 0; 0, 0, -1; 0, 1, 0]);
grad = diff()
end