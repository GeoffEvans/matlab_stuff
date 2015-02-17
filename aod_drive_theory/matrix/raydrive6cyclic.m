function sols = raydrive6cyclic()

x00 = add_ray_constraint([0;0]);
x10 = add_ray_constraint([1;0]);
x01 = add_ray_constraint([0;1]);
x11 = add_ray_constraint([1;1]);

eqs = [x00, x10, x01, x11];
eqs = simplify(eqs(:));
sols = solve(eqs,'c1','c2','c3','c4','c5','c6');
end

function [x2, k2] = prop_ray(x1, k1, w, d, D, t) % from entrance to entrance
    k2 = k1 - D * w * (dot(D,x1) - t);
    x2 = x1 + k2 * d;
end

function v = add_ray_constraint(x1)
    syms c1 c2 c3 c4 c5 c6 L f t vx vy 

    L = 0;
    C = cos((pi/3));
    S = sin((pi/3));
    D1 = [1;0];
    D2 = [C;S];
    D3 = [-C;S];
    
    k = [0;0];
    [x2,k1] = prop_ray(x1,k,c1,L,D1,t);
    [x3,k2] = prop_ray(x2,k1,c2,L,D2,t);
    [x4,k3] = prop_ray(x3,k2,c3,L,D3,t);
    [x5,k4] = prop_ray(x4,k3,c4,L,-D1,t);
    [x6,k5] = prop_ray(x5,k4,c5,L,-D2,t);
    [xF,~] = prop_ray(x6,k5,c6,f,-D3,t);
    v = diff(xF, 't')% - [vx;vy];
end