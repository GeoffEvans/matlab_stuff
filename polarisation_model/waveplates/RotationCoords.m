function [ majorAxis, Xe, Ye ] = RotationCoords( O, X, Y )

Xe = X * cos(O) - Y * sin(O);
Ye = Y * cos(O) + X * sin(O);

majorAxis = [-5*cos(O) 5*cos(O) ; -5*sin(O) 5*sin(O)];

end

