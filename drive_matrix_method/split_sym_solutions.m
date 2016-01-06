function procd = split_sym_solutions(sols)

    fields = fieldnames(sols);
    
    procd = cell(2,1);
    for m = 1:2
        temp_sym = sym(zeros(6,1));
        for k  = 1:6
            eqs = sols.(fields{k});
            temp_sym(k) = eqs(m);
        end
        procd{m} = temp_sym;
    end

end

