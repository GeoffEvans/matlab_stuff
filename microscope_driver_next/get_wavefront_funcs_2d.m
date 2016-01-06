function [w1_fun, w2_fun] = get_wavefront_funcs_2d()
    [w1_sym, w2_sym] = get_symbolic_wavefronts();
    w1_fun = matlabFunction(w1_sym);
    w2_fun = matlabFunction(w2_sym);
end 

function [w1, w2] = get_symbolic_wavefronts()
    syms C D Z0 v X0 Xpr wp1_1 wp1_2 wp1_3 wp1_4 w2_1 w2_2 w2_3 w2_4 L V F0 lambda

    x4 = w2_4 + wp1_4 - 24*D ;
    tx3 = w2_4 - wp1_4;
    x3 = w2_3 + wp1_3;
    tx1 = w2_3 - wp1_3 - C/Z0/V;
    x2 = w2_2 + wp1_2 - 1/Z0;
    tx = w2_2 - wp1_2 + v/Z0/V;
    x = w2_1 + wp1_1 + X0/Z0;
    freq = (lambda * F0 / V) - wp1_1;

    s4 = solve([x4; tx3], w2_4, wp1_4);
    s3 = solve([x3; tx1], w2_3, wp1_3);
    s2 = solve([x2; tx], wp1_2, w2_2);
    s1 = solve([x, freq], wp1_1, w2_1);

    w1 = propagate_wavefront([s1.wp1_1; s2.wp1_2; s3.wp1_3; s4.wp1_4], -L);
    w2 = [s1.w2_1; s2.w2_2; s3.w2_3; s4.w2_4];
end

function wp = propagate_wavefront(w, L)
    wp_syms = propagate_wavefront_syms();
    wp = subs(wp_syms, {'L'; 'w1'; 'w2'; 'w3'; 'w4'}, [L; num2cell(w)]);
end

function wp_0 = propagate_wavefront_syms()
    syms x w1 w2 w3 w4 L

    y = w1*x + w2*x^2/2 + w3*x^3/6 + w4*x^4/24;
    xp = x - diff(y)*L/sqrt(1+diff(y)^2);
    yp = y + L/sqrt(1+diff(y)^2) - L;

    wp = sym(zeros(4,1));
    wp(1) = diff(yp)/diff(xp);
    wp(2) = diff(wp(1))/diff(xp);
    wp(3) = diff(wp(2))/diff(xp);
    wp(4) = diff(wp(3))/diff(xp);
    %wp(5) = diff(wp(4))/diff(xp); % quartic chirp currently unsupported

    wp_0 = simplify(subs(wp, x, 0));
end


