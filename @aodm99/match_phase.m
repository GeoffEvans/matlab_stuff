function [ acAngleOpt,dAngleOpt, acInvWavelenOpt, nOrdOpt, nExt ] = match_phase( iAngleArray, acFreqArray, displayPlot )
% Finds the diffracted and acoustic beam angles for pairs of incident
% angles and acoustic frequencies.
% The pairs must be in sequence as two seperate horizontal arrays.

if nargin < 3
    displayPlot = true;
end

iAngleSize = size(iAngleArray);
acFreqSize = size(acFreqArray);
if iAngleSize(1) ~= 1 || acFreqSize(1) ~= 1 || iAngleSize(2) ~= acFreqSize(2)
    error('Inputs must be horizontal arrays of the same length.');
end

opWavelenVac = aodm99.opWavelenVac;
angularRes = 0.0005;
dif = 0.002;
[ ~, nExt, ~, ~ ] = teo2.find_n_op( abs(iAngleArray) );
[ acApproxAngle, dApproxAngle ] = CalcApproxAngles();
[acAngleOpt,dAngleOpt,acInvWavelenOpt,nOrdOpt] = CalcPreciseAngles(acApproxAngle, dApproxAngle);

if displayPlot
    Plot();
end

    function [ acApproxAngle, dApproxAngle ] = CalcApproxAngles()
        nOrdApprox = 2.2597;
        options = optimset('Algorithm','trust-region-reflective','display','off',...
            'JacobPattern',speye(length(iAngleArray)));
        acApproxAngle = fsolve(@ZeroFunction, 0*iAngleArray, options);
        invAcApproxWvln = CalcInvAcWavelen(acApproxAngle);
        nOrdSinApprox = CalcnOrdSin(acApproxAngle, invAcApproxWvln);
        dApproxAngle = asin(nOrdSinApprox / nOrdApprox);
        
        function f = ZeroFunction(acAngle)
            if isinf(acAngle)
                error('Your iAngles are too big! Have you mixed up degrees/radians?')
            end
            invAcWvln = CalcInvAcWavelen(acAngle);
            f = (CalcnOrdSin(acAngle, invAcWvln).^2 + CalcnOrdCos(acAngle, invAcWvln).^2) - nOrdApprox^2;
        end
        function nOrdSin = CalcnOrdSin(acAngle, invAcWavelen)
            nOrdSin = nExt.*sin(iAngleArray) - opWavelenVac.*invAcWavelen.*cos(acAngle);
        end
        function nOrdCos = CalcnOrdCos(acAngle, invAcWavelen)
            nOrdCos = nExt.*cos(iAngleArray) + opWavelenVac.*invAcWavelen.*sin(acAngle);
        end
    end

    function [acAngleOpt,dAngleOpt,invAcWavelenOpt,nOrdOpt] = CalcPreciseAngles(acApproxAngle, dApproxAngle)
        acAngleRanges = GetAcAngleRanges(acApproxAngle); % These should be matrices with a column of values for each input
        dAngleRanges = GetDAngleRanges(dApproxAngle);
        invAcWavelenRanges = CalcInvAcWavelen(acAngleRanges);
        [ nOrdRanges, ~, ~, ~ ] = teo2.find_n_op( abs(dAngleRanges) );
        [acAngleExtended, dAngleExtended] =	SetCrossProduct(acAngleRanges, dAngleRanges); % Get arrays containing all permutations
        [invAcWavelenExtended, nOrdExtended] =	SetCrossProduct(invAcWavelenRanges, nOrdRanges); % Get arrays containing all permutations
        
        sizeExtended = size(dAngleExtended); % All extended quantities should have same dims.
        columnOnes = ones(sizeExtended(1),1);
        % Find parameters which best satisfies (4) and (5):
        eq4 = columnOnes*(nExt.*sin(iAngleArray)) - nOrdExtended.*sin(dAngleExtended) - opWavelenVac.*invAcWavelenExtended.*cos(acAngleExtended);
        eq5 = columnOnes*(nExt.*cos(iAngleArray)) - nOrdExtended.*cos(dAngleExtended) + opWavelenVac.*invAcWavelenExtended.*sin(acAngleExtended);
        
        offset = 10 * abs(eq5) + abs(eq4);
        [~,minimiserArray] = min(offset);
        minimiser = false(size(offset));
        for m = 1:size(offset,2)
            minimiser(minimiserArray(m),m) = 1;
        end
        acAngleOpt = acAngleExtended(minimiser)';
        dAngleOpt = dAngleExtended(minimiser)';
        nOrdOpt = nOrdExtended(minimiser)';
        invAcWavelenOpt = invAcWavelenExtended(minimiser)';
        
        function acAngleRange = GetAcAngleRanges(acApproxAngle)
            rangeColumn = [(-5*dif):angularRes:0]';
            acAngleRange = rangeColumn * (1+0*acApproxAngle) + (1+0*rangeColumn) * acApproxAngle;
        end
        function dAngleRange = GetDAngleRanges(dApproxAngle)
            rangeColumn = [(-dif):angularRes:dif]';
            dAngleRange = rangeColumn * (1+0*dApproxAngle) + (1+0*rangeColumn) * dApproxAngle;
        end
    end

    function invAcWavelen = CalcInvAcWavelen(acAngle)
        theta = pi/2 - acAngle; % Sound wave goes down [110] normal to [001]
        phi = pi/4 + 0*theta;   % Fixed down [110], neglect diffraction in phi-direction
        acVel = teo2.find_v_ac_min(theta, phi);
        sizeExtended = size(acAngle); % If acAngle is scalar, no change, if acAngle is extended matrix, extend.
        columnOnes = ones(sizeExtended(1),1);
        invAcWavelen = columnOnes * acFreqArray ./ acVel;
    end

    function Plot()
        iX = nExt/opWavelenVac.*sin(iAngleArray);
        iY = nExt/opWavelenVac.*cos(iAngleArray);
        dX = nOrdOpt/opWavelenVac.*sin(dAngleOpt);
        dY = nOrdOpt/opWavelenVac.*cos(dAngleOpt);
        
        figure();
        if isscalar(iAngleArray) && isscalar(acFreqArray)
            PlotWavevectors();
        else
            PlotErrors();
        end
        
        function PlotErrors()
            xError = iX - 1.*acInvWavelenOpt.*cos(acAngleOpt) - dX;
            yError = iY + 1.*acInvWavelenOpt.*sin(acAngleOpt) - dY;
            errorSqr = xError.^2 + yError.^2;
            scatter(iAngleArray, acFreqArray, errorSqr);
            xlabel('incident angle');
            ylabel('acoustic freq');
        end
        
        function PlotWavevectors()
            polarRange = linspace(-pi/20,pi/20,2000);
            [ nPolarOrd, nPolarExt, ~, ~ ] = teo2.find_n_op( polarRange );
            plot(sin(polarRange).*nPolarOrd/opWavelenVac, cos(polarRange).*nPolarOrd/opWavelenVac);
            xlabel('[110]');
            ylabel('[001]');
            hold on;
            plot(sin(polarRange).*nPolarExt/opWavelenVac, cos(polarRange).*nPolarExt/opWavelenVac);
            line([iX*0.997, iX],[iY*0.997, iY], 'color', 'red'); % incident
            line([dX*0.997, dX],[dY*0.997, dY], 'color', 'cyan'); % diffracted
            line([iX iX - 1.*acInvWavelenOpt.*cos(acAngleOpt)],... % acoustic
                [iY iY + 1.*acInvWavelenOpt.*sin(acAngleOpt)], 'color', 'green');
            hold off;
        end
    end
end

function [ extendedX, extendedY ] = SetCrossProduct( X, Y )

% Given two matrices, X and Y, the function creates new matrices, extendedX and
% extendedY, that contain all possible permutations of the rows of X
% and Y.

% If X is of height Lx and Y is of height Ly then the extended vectors will
% be of height Lx * Ly.

sizeY = size(Y);
sizeX = size(X);
noRowsY = sizeY(1);
noRowsX = sizeX(1);
extendedX = repmat(X, sizeY(1), 1);

noRowsExtended = sizeY(1)*noRowsX;
extendedY = zeros( noRowsExtended, sizeY(2) );
for k = 1:noRowsY
    indexArray = (1:noRowsX) + noRowsX*(k-1);
    extendedY(indexArray, :) = repmat(Y(k,:),noRowsX,1);
end

end

