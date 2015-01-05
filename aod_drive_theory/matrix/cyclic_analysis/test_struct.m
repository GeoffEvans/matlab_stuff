function  test_struct(s, n)
    
    fields = fieldnames(s);   
    
    for i = 1:numel(fields)
        array = s.(fields{i});
        display(symvar(array(n)))
    end

end