wavelen = 800e-9;
l1 = 0.05;
v1 = 619;

a = parallel2();
b = subs(a, {'v1', 'l1', 'vx'}, {v1, l1, 0});
[v2,l2] = meshgrid(620:1:700,1:0.1:20);
vx = -l2 * 100;
figure()
hold on

surf(v2, l2, v1.^2/wavelen .* (v2 - vx)./(l1.*(v2-vx) + l2.*(v2-v1)),'EdgeColor','none')
freezeColors()

surf(v2, l2, -v2.^2/wavelen .* (v1 - vx)./(l2.*v1 - l2.*v2),'EdgeColor','none')
freezeColors()

surf(v2, l2, (v1.^2-v1.^3./v2)./l1./wavelen,'EdgeColor','none')
colormap hot

alpha(1);
hold off