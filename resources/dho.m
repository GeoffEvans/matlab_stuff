a = 1;
b = (1:0.8:5)';
c = 4;

L_ = (-b - sqrt(b.^2 - 4*c)) / 2;
Lp = (-b + sqrt(b.^2 - 4*c)) / 2;

A_ = a .* Lp ./ (Lp - L_);
Ap = a .* L_ ./ (L_ - Lp);

bC = 4;
L = -bC/2;

A0 = a;
A1 = -L*a;

t = 0:0.01:20;
x = Ap*(1+0*t) .* exp(Lp * t) + A_*(1+0*t) .* exp(L_ * t);
xC = (A0*(1+0*t) +  A1*t) .* exp(L * t);

plot(t,xC, 'r');
grid
hold on;
for k = 1:length(b)
    plot(t,x(k,:));
end
xlabel('time')
ylabel('position')
hold off;
