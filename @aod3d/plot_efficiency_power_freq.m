function plot_efficiency_power_freq()

acFreqMax = 100e6;
acFreqMin = 20e6;
range = 0:80;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
acPowerRange = 0.2:0.1:30;
[acFreq ,acPower ] = aod3d.SetCrossProduct( acFreqRange, acPowerRange);

tic
arrayLength = (1+0*acPower);
iThetaAir = 2.26 * pi / 180 * arrayLength;
iPhi = -3*pi/4 * arrayLength;
iInten = 1;
iPolAir = [1; 1i]/sqrt(2) * arrayLength;
acTheta = pi/2 * arrayLength;
[ ~,dThetaAir, dPhi, dInten, dPol ] = aod3d.aod_propagator( iThetaAir, iPhi, iInten, iPolAir, acFreq, acTheta, acPower );
toc

figure();
PlotEfficiencies();

    function PlotEfficiencies()
        acPowerMesh = reshape(acPower, length(acPowerRange), length(acPower)/length(acPowerRange));
        efficiencies = reshape(dInten, length(acPowerRange), length(acPower)/length(acPowerRange));
        acFreqMesh = reshape(acFreq, length(acFreq)/length(acFreqRange), length(acFreqRange));
        surf(acPowerMesh, acFreqMesh / 1e6, efficiencies, 'linestyle', 'none');
        xlabel('acoustic power / W')
        ylabel('acoustic freq / MHz')
        zlabel('efficiency')
        grid on;
        grid minor;
        axis square;
    end
end

