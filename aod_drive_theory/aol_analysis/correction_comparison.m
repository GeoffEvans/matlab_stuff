function correction_comparison()

NoCorrection6aolLong()
%GoodCorrection6aolLong()

end

function NoCorrection4aol()
correction =  0;
focalLength = 1e0;
angleSamples = 8;
aolFunction = @aol_chirps.Aod4;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function PoorCorrection4aol()
correction =  13e19;
focalLength = 1e0;
angleSamples = 8;
aolFunction = @aol_chirps.Aod4;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function GoodCorrection6aol()
correction =  10.4e19;
focalLength = 1e0;
aolFunction = @aol_chirps.Aod6pairedAS;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function NoCorrection6aol()
correction =  0;
focalLength = -1e0;
aolFunction = @aol_chirps.Aod6pairedAS;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function FocalPlane6aol()
correction =  0;
focalLength = 1e24;
aolFunction = @aol_chirps.Aod6pairedAS;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function FocalPlane4aol()
correction =  0;
focalLength = 1e24;
aolFunction = @aol_chirps.Aod4;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function NoCorrection4aolLong()
correction =  0;
focalLength = 1e1;
angleSamples = 8;
aolFunction = @aol_chirps.Aod4;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function PoorCorrection4aolLong()
correction =  2.7e19;
focalLength = 1e1;
angleSamples = 8;
aolFunction = @aol_chirps.Aod4;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function GoodCorrection6aolLong()
correction = 3.9e20;
focalLength = 1e0;
aolFunction = @aol_chirps.Aod6cyclicNewScan;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
function NoCorrection6aolLong()
correction =  0;
focalLength = 1e0;
aolFunction = @aol_chirps.Aod6cyclicNewScan;
angleSamples = 8;
aol_analysis(focalLength, correction, angleSamples, aolFunction)
end
