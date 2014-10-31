function [] = plot_fourier_series( formula )

x = 0:0.1:8;
accumulator = 0;

for n=1:10
    accumulato
    
    r = accumulator + formula(n) .* cos(pi./2.*n.*x);
end

plot(x,accumulator)

end