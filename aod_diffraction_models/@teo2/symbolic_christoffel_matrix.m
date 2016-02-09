function christoffelMatrix = symbolic_christoffel_matrix()

syms sin(t) cos(t) sin(p) cos(p)
left = [    sin(t).*cos(p)      0                   0; 
            0                   sin(t).*sin(p)      0;
            0                   0                   cos(t)];
        
right = [   0                   cos(t)              sin(t).*sin(p); 
            cos(t)              0                   sin(t).*cos(p); 
            sin(t).*sin(p)      sin(t).*cos(p)   0];
 
directionMatrix = [left, right];
christoffelMatrix = directionMatrix * teo2.ElasticStiffness() * transpose(directionMatrix);
end