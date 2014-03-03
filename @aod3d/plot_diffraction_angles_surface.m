function plot_diffraction_angles_surface()

acFreqMax = 60e6;
acFreqMin = 20e6;
range = 0:80;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
iThetaAirDegRange = 2.1:0.01:2.4;
iThetaAirRange = iThetaAirDegRange / 180 * pi;
[acFreq ,iThetaAir ] = aod3d.SetCrossProduct( acFreqRange, iThetaAirRange);
tic
arrayLength = (1+0*iThetaAir);
iPhi = -3*pi/4 * arrayLength;
iInten = 1;
iPolAir = [1; 1i]/sqrt(2) * arrayLength;
acTheta = pi/2 * arrayLength;
acPower = 1.8;
[ ~,dThetaAir, dPhi, dInten, dPol ] = aod3d.aod_propagator( iThetaAir, iPhi, iInten, iPolAir, acFreq, acTheta, acPower );
toc

figure();
absTheta = IncludeNegativeThetas(dThetaAir,dPhi)-iThetaAir;
PlotDiffractionAngles(absTheta, 'absolute');

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

















