function plot_efficiency_surface_2aod()

acFreqMax = 125e6;
acFreqMin = 5e6;
range = 0:30;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
iThetaAirDegRange = 1.6:0.04:4.2;
iThetaAirRange = iThetaAirDegRange / 180 * pi;
[acFreq ,iThetaAir ] = aod3dFirst.SetCrossProduct( acFreqRange, iThetaAirRange);
tic
iPhi = -3*pi/4 * 1+0*iThetaAir;
[dThetaAir, dPhi, dIntensity, dPolarisation] = aod3dFirst.aodPropagator(iThetaAir, iPhi, acFreq, 1, [1; 1i] * (1/sqrt(2)+0*iThetaAir));
polEff = PutThroughQwpLpQwp();
[dThetaAir, dPhi, dIntensity, dPolarisation] = aod3dFirst.aodPropagator(iThetaAir, iPhi, acFreq, dIntensity, dPolarisation);
toc

PlotOverallEfficiency();

    function polEff = PutThroughQwpLpQwp()
        % qwp: quarter wave plate
        qwpPol = dPolarisation;
        qwpPol(2,:) = dPolarisation(2,:)*1i;
        linearPol = dPolarisation*0 + 1/sqrt(2);
        polComponent = dot(qwpPol, linearPol);
        polEff = polComponent .* conj(polComponent);
        dIntensity = dIntensity .* polEff;
        dPolarisation = qwpPol;
        dPolarisation(2,:) = qwpPol(2,:)*1i;
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















