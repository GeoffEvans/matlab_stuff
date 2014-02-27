function [ velocities, displacementVecs ] = find_v_ac_sorted( thetaMesh, phiMesh )
% Returns sorted cell arrays containing three modes' velocities and
% displacement vectors.

[eVals,dispVecs] = teo2.find_v_ac_inner( thetaMesh, phiMesh );
[displacementVecs, eVals] = cellfun(@GetSortedEigenvalues, dispVecs, eVals, 'UniformOutput', false);
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

    function [eVecsSorted, eValsSorted] = GetSortedEigenvalues(eVecs, eVals)
        [eValsSorted, indicesOfSortedArray]=sort(eVals);
        eVecsSorted = eVecs(:,indicesOfSortedArray); % arrange the columns in this order
    end
end