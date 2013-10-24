% Calculate the plane wave interaction efficiency using (10) from Maak 1999

acFreqMax = 120e6;
acFreqMin = 5e6;

range = 1:200;
modelVals = range;
for k = [0.8,0.9,1,1.1,1.5,2,3]
    for j = range;
        freq = acFreqMin + j/max(range)*acFreqMax;
        modelVals(j) = maak_eff(freq, k);
    end

    freqRange = range/max(range)*acFreqMax;
    plot(freqRange, modelVals, 'color', rand(1,3));
    hold on;
    grid on;
    grid minor;
end