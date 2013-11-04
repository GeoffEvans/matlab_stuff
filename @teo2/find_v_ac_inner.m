function [ eVals, displacementVectors ] = find_v_ac_inner( thetaMesh, phiMesh )
% Returns unsorted cell arrays containing three modes' eigenvalues and
% eigenvectors.

% Solve simple eigenvalue equation (1.107) Xu & Stroud:
christoffelMatrices = arrayfun(@ChristoffelMatrix, thetaMesh, phiMesh, 'UniformOutput', false);
[displacementVectors, eVals] = cellfun(@GetEigenvaluesAsArray, christoffelMatrices, 'UniformOutput', false);

    function f = ChristoffelMatrix(t,p) % Ugly but optimised
        f = [ 26500000000*cos(t)^2 + 55700000000*cos(p)^2*sin(t)^2 + 65900000000*sin(p)^2*sin(t)^2, 117100000000*cos(p)*sin(p)*sin(t)^2, 48300000000*cos(p)*cos(t)*sin(t);
            117100000000*cos(p)*sin(p)*sin(t)^2, 26500000000*cos(t)^2 + 65900000000*cos(p)^2*sin(t)^2 + 55700000000*sin(p)^2*sin(t)^2, 48300000000*cos(t)*sin(p)*sin(t);
            48300000000*cos(p)*cos(t)*sin(t), 48300000000*cos(t)*sin(p)*sin(t), 105800000000*cos(t)^2 + 26500000000*cos(p)^2*sin(t)^2 + 26500000000*sin(p)^2*sin(t)^2];
    end

    function [eVecs, eVals] = GetEigenvaluesAsArray(matrix)
        [eVecs, eValsMatrix] = eig(matrix);
        eVals = diag(eValsMatrix);
    end
end

