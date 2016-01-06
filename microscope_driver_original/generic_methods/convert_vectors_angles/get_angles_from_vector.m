function [ modulus, theta, phi ] = get_angles_from_vector( vectorMatrix )

modulus = mag( vectorMatrix );
phi = atan2( vectorMatrix(2,:,:), vectorMatrix(1,:) );
theta = acos( vectorMatrix(3,:,:) ./ modulus ) ; % Assumes WLOG theta > 0

end
