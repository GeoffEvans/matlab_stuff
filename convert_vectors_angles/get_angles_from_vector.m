function [ modulus, theta, phi ] = get_angles_from_vector( vectorMatrix )

modulus = GetModulusColumns( vectorMatrix );
phi = atan2( vectorMatrix(2,:), vectorMatrix(1,:) );
theta = acos( vectorMatrix(3,:) ./ modulus ) ; % Assumes WLOG theta > 0

end

function modulus = GetModulusColumns( matrix )
    modulusSqr = sum( matrix .^ 2 );
    modulus = sqrt( modulusSqr );
end
