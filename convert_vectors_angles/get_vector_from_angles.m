function [ vector ] = get_vector_from_angles( modulus, theta, phi )

x = modulus .* sin(theta) .* cos(phi);
y = modulus .* sin(theta) .* sin(phi);
z = modulus .* cos(theta);

vector = [x;y;z];

end
