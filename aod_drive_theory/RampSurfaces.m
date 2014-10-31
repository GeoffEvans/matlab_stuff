wavelen = 800e-9;
l1 = 0.05;
v1 = 619;

a = parallel2();
b = subs(a, {'v1', 'l1', 'vx'}, {v1, l1, 0});
[vx,l2] = meshgrid(-5:1:5,1:0.1:20);
figure()
hold on

colormap cool
surf(vx * v1, l2, abs(v1.^2/wavelen .* (vx + 1)./(l1 + 2*l2 + l1.*vx)),'EdgeColor','none')
freezeColors()

surf(vx * v1, l2, abs(v1.^2/wavelen .* -(vx - 1)./(2*l2)),'EdgeColor','none')
colormap summer

xlabel('Vx')
ylabel('L2')
hold off