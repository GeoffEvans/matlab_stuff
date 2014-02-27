function [ p ] = PlotEllipse( X, Y )

p = plot(X,Y, 'color', 'green');
xName = inputname(1);
yName = inputname(2);
set(p,'XDataSource',xName);
set(p,'YDataSource',yName);
grid on;
grid minor;
axis([-5 5 -5 5]);
axis square;

end

