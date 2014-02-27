constT2 = @(b,x,t2) b.*x + t2;
constX2 = @(b,x,x2) (x - x2) ./ b;
xMax = 10;
xMin = -2;
X = xMin:xMax;
b = 0.3;
plot(X, constT2(b,X,0),'m','LineWidth',2);
hold on
plot(X, constX2(b,X,0),'m','LineWidth',2);
axis([xMin xMax xMin xMax])
axis square
plot([xMin xMax],[0 0],'LineWidth',2,'color','black');
plot([0 0], [xMin xMax],'LineWidth',2,'color','black');

for k = (xMin-5):(xMax+5)
    plot(X, constT2(b,X,k), ':m');
    plot(X, constX2(b,X,k), ':m');
    plot([xMin xMax],[k k], ':k');
    plot([k k], [xMin xMax],':k');
end
hold off;