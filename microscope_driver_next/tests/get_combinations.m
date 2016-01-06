function [ combs ] = get_combinations(vec)
[p,q] = meshgrid(vec, vec);
combs = [p(:)'; q(:)'];
end

