function saw_fourier_series()

a = 2;
L = 2;
X = -5:0.02:5;

    function t = term(x,n)
        coefficient = a .* (-1).^(n'+1) ./ (pi .* n') * (x*0+1); 
        phase = sin(2*pi/L .* n' * x);
        t = coefficient .* phase;
    end

    function s = series(x)
        s = sum(term(x,1:100));
    end

plot(X,series(X));
hold on;
plot(X,sawtooth((X-1)*pi), 'color', 'red');
hold off;
end

