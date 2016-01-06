function plot_phase_surface()

acFreqMax = 90e6;
acFreqMin = 10e6;
range = 0:80;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
iThetaAirDegRange = 0.8:0.04:2;
iThetaAirRange = iThetaAirDegRange / 180 * pi;
[acFreq ,iThetaAir ] = aod3d.SetCrossProduct( acFreqRange, iThetaAirRange);

tic
arrayLength = (1+0*iThetaAir);
iPhi = -3*pi/4 * arrayLength;
acTheta = pi/2 * arrayLength;
[~, ~, phase] = aod3d.match_phase(iThetaAir, iPhi , acFreq, acTheta );
toc

PlotOverallEfficiency();

    function PlotOverallEfficiency()
        figure();
        iThetaAirDegMesh = reshape(iThetaAir * 180/pi, length(iThetaAirRange), length(iThetaAir)/length(iThetaAirRange));
        efficiencies = reshape(abs(phase), length(iThetaAirRange), length(iThetaAir)/length(iThetaAirRange));
        acFreqMesh = reshape(acFreq, length(acFreq)/length(acFreqRange), length(acFreqRange));
        surf(iThetaAirDegMesh, acFreqMesh / 1e6, efficiencies, 'linestyle', 'none');
        xlabel('incidence angle air / degrees')
        ylabel('acoustic freq / MHz')
        zlabel('efficiency')
        grid on;
        grid minor;
        axis square;
    end
end













