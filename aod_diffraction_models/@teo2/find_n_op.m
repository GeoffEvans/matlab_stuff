function [ nOrd, nExt, pOrd, pExt ] = find_n_op( thetaRange )

% Properties of paratellurite
gammaRange = teo2.activityVector * cos(thetaRange).^2; % Xu&St (1.74)
transverseImpermEigs = FindTransverseImpermittivityEigenvals(thetaRange);

[nOrd, nExt] = FindRefractiveIndex(gammaRange,transverseImpermEigs);
[pOrd,pExt] = FindPolarisations(gammaRange,transverseImpermEigs);

% Note the polarisation DyDx does not refer to the original XY axes of the
% incident wave but to the principal axes of the transverse impermeability.

    function transverseImpermEigs = FindTransverseImpermittivityEigenvals(thetaRange)
        % The coordinates of the crystal properties are the natural coords of the crystal.
        % Rotate the YZ-axes of this coord system so that the new Z-axis orientation is given by incident optic wave.
        transverseImpermEigs = cell(1,length(thetaRange));
        for q = 1:length(thetaRange)
            transverseImpermEigs{q} = TransverseImpermEigenvals(thetaRange(q));
        end
    end
    function [nOrd, nExt] = FindRefractiveIndex(gammaRange,transverseImpermEigs)
        % Use equation (1.62) Xu & Stroud to find values for n for this new coord system.        
        Disc = @(eigVals, gamma) (eigVals(2) - eigVals(1))^2 + 4*gamma^2;
        CalculateN1 = @(eigVals, discs) (0.5 * ( eigVals(1) + eigVals(2) - sqrt(discs)))^-0.5;
        CalculateN2 = @(eigVals, discs) (0.5 * ( eigVals(1) + eigVals(2) + sqrt(discs)))^-0.5;
        
        discrims = cellfun(Disc, transverseImpermEigs, num2cell(gammaRange), 'UniformOutput', false);
        nExt = cellfun(CalculateN1, transverseImpermEigs, discrims);
        nOrd = cellfun(CalculateN2, transverseImpermEigs, discrims);
    end
    function [pOrd,pExt] = FindPolarisations(gammaRange,transverseImpermEigs)
        transverseImpermEigs = cell2mat(transverseImpermEigs);
        pOrd = 1i./(gammaRange) .* (transverseImpermEigs(1,:) - nOrd.^-2);
        pExt = 1i./(gammaRange) .* (transverseImpermEigs(1,:) - nExt.^-2);
    end
end

function transverseEigenVals = TransverseImpermEigenvals(theta)
    yzRotationMat = [1 0 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];
    transverseImpermeability = yzRotationMat * teo2.RelativeImpermeability() * (yzRotationMat'); % (1.59) diagonalise XY components of transverseImpermeability (not eigenvalues of n).
    transverseEigenVals = diag(transverseImpermeability(1:2,1:2));
end

