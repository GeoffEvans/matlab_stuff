function correction_comparison()

NoCorrection6aol()
GoodCorrection6aol()

end

function NoCorrection4aol()
correction =  0;
focalLength = 1e1;
angleSamples = 8;
aolFunction = @aol_chirps.Aod4;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function PoorCorrection4aol()
correction =  2.7e19;
focalLength = 1e1;
angleSamples = 8;
aolFunction = @aol_chirps.Aod4;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function NoCorrection6aol()
correction =  0;
focalLength = 1e1;
aolFunction = @aol_chirps.Aod6pairedAS;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function GoodCorrection6aol()
correction =  2.3e19;
focalLength = 1e1;
aolFunction = @aol_chirps.Aod6pairedAS;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end