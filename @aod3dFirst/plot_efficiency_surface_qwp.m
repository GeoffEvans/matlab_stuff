function plot_efficiency_surface_qwp()

acFreqMax = 60e6;
acFreqMin = 10e6;
range = 0:20;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
iThetaAirDegRange = 1.5:0.2:4;
iThetaAirRange = iThetaAirDegRange / 180 * pi;
[acFreq ,iThetaAir ] = aod3dFirst.SetCrossProduct( acFreqRange, iThetaAirRange);
tic
phi = -3*pi/4;
[dThetaAir, dPhi, dIntensity, dPolarisation] = aod3dFirst.aod_propagator(iThetaAir, phi * 1+0*iThetaAir, acFreq, 1, [1; 1i] * (1/sqrt(2)+0*iThetaAir));
%polEff = PutThroughQwp();
toc

PlotOverallEfficiency();

    function polEff = PutThroughQwp()
        % qwp: quarter wave plate
        qwpPol = dPolarisation;
        qwpPol(2,:) = dPolarisation(2,:)*1i;
        linearPol = dPolarisation*0 + 1/sqrt(2);
        polComponent = dot(qwpPol, linearPol);
        polEff = polComponent .* conj(polComponent);
        dIntensity = dIntensity .* polEff;
    end

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
        grid minor;
        axis square;
    end
end













