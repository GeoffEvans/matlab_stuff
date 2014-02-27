function [ velocities, displacementVectors ] = find_v_ac( thetaMesh, phiMesh )
% Returns unsorted cell arrays containing three modes' velocities and
% displacement vectors.

[eVals, displacementVectors] = teo2.find_v_ac_inner( thetaMesh, phiMesh );
velocities = GetVelocitiesFromEvals(eVals);

    function velocities = GetVelocitiesFromEvals(eValsCells)
        velocitiesMatrix = GetVelFromEigs(cell2mat(eValsCells));
        velMatSize = size(velocitiesMatrix);
        velocitiesMatrixReshaped = reshape(velocitiesMatrix,3,velMatSize(1)/3 * velMatSize(2));
        velocitiesCellUnshaped = num2cell(velocitiesMatrixReshaped,1);
        velocities = reshape(velocitiesCellUnshaped, velMatSize(1)/3, velMatSize(2));
        
        function f = GetVelFromEigs(eigs)
            f = sqrt(eigs/teo2.density);
        end
    end
end
