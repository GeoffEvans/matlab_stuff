function shifted_coeffs = poly_shift(coeffs, shift)
    % expects list of vectors [(a;b;c;...),...]
    
    shifted_coeffs = zeros(size(coeffs));
    num_coeffs = size(coeffs, 1);
    num_polys = size(coeffs, 2);
    if numel(shift) == 1
        shift = repmat(shift, 1, num_polys);
    end
    
    for n = 1:num_coeffs
        binom = [diag(fliplr(pascal(n))); zeros(num_coeffs - n, 1)];
        shift_powers = [repmat(shift, n, 1) .^ repmat(fliplr(0:(n-1))', 1, num_polys); zeros(num_coeffs - n, num_polys)];
        n_coeffs = coeffs(n,:);

        shifted_coeffs = shifted_coeffs + (binom * n_coeffs) .* shift_powers;
    end
end
