function [ nOrd, nExt, pOrd, pExt ] = find_n_op( thetaRange )

% Properties of paratellurite
gammaRange = teo2.activityVector * cos(thetaRange).^2; % Xu&St (1.74)
transverseEigenVals = FindTransversePermittivityEigenvals();
[nOrd, nExt] = FindN();
[pOrd,pExt] = FindPolarisations();
% Note the polarisation DyDx does not refer to the original XY axes of the
% incident wave but to the principal axes of the transverse impermeability.

    function transverseEigenVals = FindTransversePermittivityEigenvals()
        % The coordinates of the crystal properties are the natural coords of the crystal.
        % Rotate the YZ-axes of this coord system so that the new Z-axis orientation is given by incident optic wave.
        % This assumes the phi is pi/2. 
        transverseEigenVals = cell(1,length(thetaRange));
        for q = 1:length(thetaRange)
            transverseEigenVals{q} = RefIndEigenvals(thetaRange(q));
        end
    end
    function [nOrd, nExt] = FindN()
        % Use equation (1.62) Xu & Stroud to find values for n for this new coord system.        
        Disc = @(eigVals, gamma) (eigVals(2) - eigVals(1))^2 + 4*gamma^2;
        CalculateN1 = @(eigVals, discs) (0.5 * ( eigVals(1) + eigVals(2) - sqrt(discs)))^-0.5;
        CalculateN2 = @(eigVals, discs) (0.5 * ( eigVals(1) + eigVals(2) + sqrt(discs)))^-0.5;
        discrims = cellfun(Disc, transverseEigenVals, num2cell(gammaRange), 'UniformOutput', false);
        nExt = cellfun(CalculateN1, transverseEigenVals, discrims);
        nOrd = cellfun(CalculateN2, transverseEigenVals, discrims);
    end
    function [pOrd,pExt] = FindPolarisations()
        pOrd = cellfun(@CalculateDyDx, num2cell(nOrd), transverseEigenVals, num2cell(gammaRange));
        pExt = cellfun(@CalculateDyDx, num2cell(nExt), transverseEigenVals, num2cell(gammaRange));
        
        function dydx = CalculateDyDx(n, transverseEigVals, gamma)
            % Taking values for n, calculate the polarisation, Dy/Dx using (1.60)
            % These will be complex (model uses complex waveforms with phase {wt - k.r})
            dydx = 1i./(gamma) .* (transverseEigVals(1) - (n)^-2);
        end
    end
end

function transverseEigenVals = RefIndEigenvals(theta)
    yzRotationMat = [1 0 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];
    transverseImpermeability = yzRotationMat * teo2.RelativeImpermeability() * (yzRotationMat'); % (1.59) diagonalise XY components of transverseImpermeability (not eigenvalues of n).
    transverseEigenVals = diag(transverseImpermeability(1:2,1:2));
end

