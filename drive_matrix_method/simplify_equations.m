function m = simplify_equations(eqs)
    m = expand(eqs);
    m = subs(m, 'L^5', 0);
    m = subs(m, 'L^4', 0);
    m = subs(m, 'L^3', 0);
    m = subs(m, 'L^2', 0);
    
    pwr = 2;
    coms = nchoosek(1:6, pwr);
    for n = 1:size(coms, 1)
        str = ['c', num2str(coms(n,1))];
        for nn = 2:pwr
            str = [str, '*c', num2str(coms(n,nn))];
        end
        m = subs(m, str, 0);
    end
    m = simplify(m);
end