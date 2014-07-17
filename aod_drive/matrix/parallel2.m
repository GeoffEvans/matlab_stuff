function sols = parallel2()

syms c1 c2 l1 l2 v1 v2 vx
I2 = [1 0; 0 1];

P1 = [1 l1; 0 1];
P2 = [1 l2; 0 1];

Q1 = [1 0; -c1 1];
Q2 = [1 0; -c2 1];

M = P2*Q2*P1*Q1;
%  M(1,1) == 0

D12 = P2*(Q2*P1*c1*v1+I2*c2*v2);

D = D12(1,2) - vx;

eqs = [M(1,1) D];
eqs = eqs(:);
s = solve(eqs,'c1','c2');
sols = [s.c1; s.c2];

end

