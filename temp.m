n = 1000;
u0 = 3.5;
m = linspace(0.5,0.9,n);
options = optimset('display','off');

tic
x1 = arrayfun(@(m)fsolve(@(u)ellipj(u,m),u0,options),m);
t1 = toc;

tic
x2 = zeros(1,n);
for i = 1:n
    x2(i) = fsolve(@(u) ellipj(u,m(i)),u0,options);
end
t2 = toc;

u0 = u0*ones(1,n);
tic
options = optimset('Algorithm','trust-region-reflective','display','off','largescale','on',...
    'JacobPattern',speye(n));
x3 = fsolve(@(u,m) ellipj(u,m),u0,options,m);
t3 = toc;

[x1' x2' x3']
[t1 t2 t3]