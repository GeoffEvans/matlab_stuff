p = perms([2,3,4,5,6]);
order = [ones(48,1),  p(p(:,1) < 4, :)];

sols = cell(48, 1);
for n = 1:48
    sols{n} = matrixdrive6_136254(order(n,:));
end