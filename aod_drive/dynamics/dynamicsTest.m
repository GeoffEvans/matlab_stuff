function [ sum ] = dynamicsTest()

syms q T x y

a = [ -q/4, q/4, -q, q, -q, q];

n = 1;
C1 = sym(cos((n-1)*pi/3));
S1 = sym(sin((n-1)*pi/3));
t1 = T - x*C1 - y*S1;
C2 = sym(cos((n-1)*pi/3 + pi));
S2 = sym(sin((n-1)*pi/3 + pi));
t2 = T - x*C2 - y*S2;
n = 3;
C3 = sym(cos((n-1)*pi/3));
S3 = sym(sin((n-1)*pi/3));
t3 = T - x*C3 - y*S3;
C4 = sym(cos((n-1)*pi/3 + pi));
S4 = sym(sin((n-1)*pi/3 + pi));
t4 = T - x*C4 - y*S4;
n = 5;
C5 = sym(cos((n-1)*pi/3));
S5 = sym(sin((n-1)*pi/3));
t5 = T - x*C5 - y*S5;
C6 = sym(cos((n-1)*pi/3 + pi));
S6 = sym(sin((n-1)*pi/3 + pi));
t6 = T - x*C6 - y*S6;

t = [t1;t2;t3;t4;t5;t6];
t3 = t.^3;

sum = latex(expand(a*t3));

end

