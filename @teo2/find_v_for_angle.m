function [ velocities, eVecs ] = find_v_for_angle( thetaMesh, phiMesh )

% Properties of paratellurite
density = 5990; %kg/m3
elasticStiffness = [[5.57 5.12 2.18; 5.12 5.57 2.18; 2.18 2.18 10.58],...
    zeros(3); zeros(3) diag([2.65 2.65 6.59])] * 1e10;

% Solve simple eigenvalue equation (1.107) Xu & Stroud:
christoffelMatrices = arrayfun(@ChristoffelMatrix, thetaMesh, phiMesh, 'UniformOutput', false);
[eVecs, eVals] = cellfun(@GetSortedEigenvalues, christoffelMatrices, 'UniformOutput', false);
velocities = cellfun(@GetVelFromEigs, eVals, 'UniformOutput', false);

    function f = ChristoffelMatrix(t,p)
        directionMat = DirectionMatrix(t,p);
        f = directionMat * elasticStiffness * directionMat';
    end

    function f = GetVelFromEigs(eigs)
        f = sqrt(diag(eigs)/density);
    end
end

function [eVecs, eVals] = GetSortedEigenvalues(matrix)
    [eVecs, eVals] = eig(matrix);
    [eVecs, eVals] = sortem(eVecs,eVals);
end

function f = SymmetricMatrix(theta, phi)
f = [0 cos(theta) sin(theta).*sin(phi); cos(theta) 0 sin(theta).*cos(phi); sin(theta).*sin(phi) sin(theta).*cos(phi) 0];
end

function f = DirectionMatrix(t,p)
f = [diag([sin(t).*cos(p) sin(t).*sin(p) cos(t)]), SymmetricMatrix(t,p)];
end