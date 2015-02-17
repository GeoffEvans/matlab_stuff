function [ velocities ] = find_v_ac_min( thetaMesh, phiMesh )
% Returns the velocity of the slow mode.

eVals = teo2.find_v_ac_inner( thetaMesh, phiMesh );
velocities = GetMinVelocityFromEval(eVals);

    function velocities = GetMinVelocityFromEval(eValsCells)
        eigsMatrix = cell2mat(eValsCells);
        eigMatSize = size(eigsMatrix);
        eigsMatrixReshaped = reshape(eigsMatrix,3,eigMatSize(1)/3 * eigMatSize(2));
        minEigs = min(eigsMatrixReshaped);
        minVels = GetVelFromEigs(minEigs);
        velocities = reshape(minVels, eigMatSize(1)/3, eigMatSize(2));
        
        function f = GetVelFromEigs(eigs)
            f = sqrt(eigs/teo2.density);
        end
    end
end

