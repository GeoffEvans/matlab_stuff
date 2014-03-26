function plot_efficiency_surface()

acFreqMax = 70e6;
acFreqMin = 10e6;
range = 0:40;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
iThetaAirDegRange = 1.5:0.04:4;
iThetaAirRange = iThetaAirDegRange / 180 * pi;
[acFreq ,iThetaAir ] = aod3d.SetCrossProduct( acFreqRange, iThetaAirRange);

tic
arrayLength = (1+0*iThetaAir);
iPhi = -3*pi/4 * arrayLength;
iInten = 1;
iPolAir = [1; 1i]/sqrt(2) * arrayLength;
acTheta = pi/2 * arrayLength;
acPower = 1.8;
[~,~,~,dIntensity, dPol] = aod3d.aod_propagator(iThetaAir, iPhi , iInten, iPolAir, acFreq, acTheta, acPower );
toc

PlotOverallEfficiency();

    function PlotOverallEfficiency()
        figure();
        iThetaAirDegMesh = reshape(iThetaAir * 180/pi, length(iThetaAirRange), length(iThetaAir)/length(iThetaAirRange));
        efficiencies = reshape(dIntensity, length(iThetaAirRange), length(iThetaAir)/length(iThetaAirRange));
        acFreqMesh = reshape(acFreq, length(acFreq)/length(acFreqRange), length(acFreqRange));
        surf(iThetaAirDegMesh, acFreqMesh / 1e6, efficiencies, 'linestyle', 'none');
        xlabel('incidence angle air / degrees')
        ylabel('acoustic freq / MHz')
        zlabel('efficiency')
        grid on;
        colorbar;
        axis square;
    end
end













