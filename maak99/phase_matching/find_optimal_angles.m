function [ acAngleOpt,dAngleOpt, acWavelenOpt, nOrdOpt, nExt ] = find_optimal_angles( iAngle, angularRes, acFreq, opWavelenVac )

[ ~, nExt, ~, ~, ~ ] = find_n_for_angle( abs(iAngle) );
[ acApproxAngle, dApproxAngle ] = CalcApproxAngles();
[acAngleRange,dAngleRange] = GetSearchRanges();
[acAngleMesh, dAngleMesh] = meshgrid(acAngleRange, dAngleRange);
acWavelenMesh = GetAcWavelenMesh();
nOrdMesh = GetNordMesh();

eq4 = nExt.*sin(iAngle) - nOrdMesh.*sin(dAngleMesh) - opWavelenVac./acWavelenMesh.*cos(acAngleMesh);
eq5 = nExt.*cos(iAngle) - nOrdMesh.*cos(dAngleMesh) + opWavelenVac./acWavelenMesh.*sin(acAngleMesh);

% Find parameters which best satisfies (4) and (5):
offset = 10 * abs(eq5) + abs(eq4);
minimiser = (offset == min(offset(:)));
acAngleOpt = acAngleMesh(minimiser);
dAngleOpt = dAngleMesh(minimiser);
nOrdOpt = nOrdMesh(minimiser);
acWavelenOpt = acWavelenMesh(minimiser);

%PlotWavevectors();

    function [ acApproxAngle, dApproxAngle ] = CalcApproxAngles()
        nOrdApprox = 2.2597;
        acApproxAngle = fzero(@ZeroFunction,0);
        acApproxWvln = CalcAcWavelen(acApproxAngle);
        nOrdSinApprox = CalcnOrdSin(acApproxAngle, acApproxWvln);
        dApproxAngle = asin(nOrdSinApprox / nOrdApprox);
        
        function f = ZeroFunction(acAngle)
            acWvln = CalcAcWavelen(acAngle);
            f = (CalcnOrdSin(acAngle, acWvln).^2 + CalcnOrdCos(acAngle, acWvln).^2) - nOrdApprox^2;
        end
        
        function nOrdSin = CalcnOrdSin(acAngle, acWavelen)
            nOrdSin = nExt.*sin(iAngle) - opWavelenVac./acWavelen.*cos(acAngle);
        end
        
        function nOrdCos = CalcnOrdCos(acAngle, acWavelen)
            nOrdCos = nExt.*cos(iAngle) + opWavelenVac./acWavelen.*sin(acAngle);
        end
    end

    function acWavelen = CalcAcWavelen(acAngle)
        theta = pi/2 - acAngle; % Sound wave goes down [110] normal to [001]
        phi = acAngle*0 + pi/4;   % Fixed down [110], neglect diffraction in phi-direction
        acVel = cellfun(@min, find_v_for_angle(theta, phi));
        acWavelen = acVel ./ acFreq;
    end

    function [acAngleRange,dAngleRange] = GetSearchRanges()
        dif = 0.002;
        acAngMin = min(acApproxAngle) - 2 * dif;
        acAngMax = max(acApproxAngle) + dif;
        dAngMin = min(dApproxAngle) - dif;
        dAngMax = max(dApproxAngle) + dif;
        acAngleRange = (acAngMin:angularRes:acAngMax);
        dAngleRange = (dAngMin:angularRes:dAngMax);
    end

    function nOrdMesh = GetNordMesh()
        [ nOrdRange, ~, ~, ~, ~ ] = find_n_for_angle( abs(dAngleRange) );
        nOrdMesh = repmat(nOrdRange,length(acAngleRange),1)';
    end

    function acWavelenMesh = GetAcWavelenMesh()
        acWavelenRange = CalcAcWavelen(acAngleRange);
        acWavelenMesh = repmat(acWavelenRange,length(dAngleRange),1);
    end

    function PlotWavevectors()
        subplot(1,2,1); % Plot the intersections of the three planes to solve
        surf(acAngleMesh, dAngleMesh, dAngleMesh*0, 'EdgeColor','none','FaceColor',[0 0 0]);
        hold on;
        surf(acAngleMesh, dAngleMesh, eq4, 'EdgeColor','none','FaceColor',[1 0 0]); % red
        surf(acAngleMesh, dAngleMesh, eq5, 'EdgeColor','none','FaceColor',[0 0 1]); % blue
        xlabel('acoustic angle');
        ylabel('diffraction angle');
        hold off;
        
        subplot(1,2,2); % Plot the wavevector diagram
        polarRange = linspace(-pi/20,pi/20,2000);
        [ nPolarOrd, nPolarExt, ~, ~, ~ ] = find_n_for_angle( polarRange );
        plot(sin(polarRange).*nPolarOrd/opWavelenVac, cos(polarRange).*nPolarOrd/opWavelenVac);
        xlabel('[110]');
        ylabel('[001]');
        hold on;
        plot(sin(polarRange).*nPolarExt/opWavelenVac, cos(polarRange).*nPolarExt/opWavelenVac);
        iX = nExt/opWavelenVac*sin(iAngle);
        iY = nExt/opWavelenVac*cos(iAngle);
        dX = nOrdOpt/opWavelenVac*sin(dAngleOpt);
        dY = nOrdOpt/opWavelenVac*cos(dAngleOpt);
        line([iX*0.997, iX],[iY*0.997, iY], 'color', 'red');
        line([dX*0.997, dX],[dY*0.997, dY], 'color', 'red');      
        line([iX iX - 1/acWavelenOpt*cos(acAngleOpt)],[iY iY + 1/acWavelenOpt*sin(acAngleOpt)], 'color', 'green');
        hold off;
    end

end

