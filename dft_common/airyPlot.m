spacing = 0.1;
p = 0:spacing:1000;
r = 0:0.0001:0.03;
[mr,mp] = meshgrid(r,p);
grand = besselj(0,mp.*mr).*mp*spacing;
a = sum(grand,1);
figure()
plot(r,a)
t = -pi:0.001*pi:pi;
[tt,tr] = meshgrid(t,r);
[~,ta] = meshgrid(t,a);
s = pcolor(tr.*cos(tt), tr.*sin(tt), real(ta.^0.5));
set(s,'LineStyle', 'none')