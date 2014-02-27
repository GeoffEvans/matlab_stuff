function plot_diffraction_angles_surface()

acFreqMax = 50e6;
acFreqMin = 25e6;
range = 0:50;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
iThetaAirDegRange = 2.1:0.01:2.4;
iThetaAirRange = iThetaAirDegRange / 180 * pi;
[acFreq ,iThetaAir ] = aod3dFirst.SetCrossProduct( acFreqRange, iThetaAirRange);
tic
iPhi = -3*pi/4;
[dThetaAir, dPhi, dIntensity, dPolarisation] = aod3dFirst.aod_propagator(iThetaAir, iPhi, acFreq, 1, [1; 1i] * (1/sqrt(2)+0*iThetaAir));
toc

figure();
%subplot(1,2,1);
absTheta = IncludeNegativeThetas(dThetaAir,dPhi)-iThetaAir;
PlotDiffractionAngles(absTheta, 'absolute');
%subplot(1,2,2);
%PlotDiffractionAngles(dPhi, 'phi');

    function PlotDiffractionAngles(angles, nameAngle)
        iThetaAirDegMesh = reshape(iThetaAir * 180/pi, length(iThetaAirRange), length(iThetaAir)/length(iThetaAirRange));
        anglesMesh = reshape(angles * 180/pi, length(iThetaAirRange), length(iThetaAir)/length(iThetaAirRange));
        acFreqMesh = reshape(acFreq, length(acFreq)/length(acFreqRange), length(acFreqRange));
        surf(iThetaAirDegMesh, acFreqMesh / 1e6, anglesMesh, 'linestyle', 'none');
        xlabel('incidence angle air / degrees')
        ylabel('acoustic freq / MHz')
        zlabel(nameAngle)
        grid on;
        grid minor;
        axis square;
    end
    function thetaOut = IncludeNegativeThetas(theta,phi)
        phiNeg = (abs(phi-iPhi)<1e-3).*2-1;
        thetaOut = theta .* phiNeg;
    end
end

















