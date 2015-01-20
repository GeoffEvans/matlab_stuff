order1 = [1,2,3,4,5,6];
order2 = [1,3,2,4,6,5];
order3 = [1,2,4,6,5,3];
order4 = [1,2,5,6,3,4]; % no simple sol
order5 = [1,2,5,3,6,4];
order6 = [1,2,3,5,6,4];

order = {order1, order2, order3, order5};

N = numel(order);
sols = cell(N,1);
for n = 1:N
    sols{n} = matrixdrive6_136254(order{n});
end