function [ freqDeflect ] = freqDeflectEqu( scanVar,systemVar,apparentDepth,d4AppInUse,x,y )

%freqDeflect is computed for a general case, x and y are normalised
%coorinates, d4AppInUse(d4Appd,4AppMin,d4AppMid) changes depending on the imagingMode.

freqDeflect(1)= -systemVar.diffractionMode  *2 * (apparentDepth.d3App + d4AppInUse) * ...
    apparentDepth.d4App / ((d4AppInUse + apparentDepth.d3App) * ...
    (apparentDepth.d1App + apparentDepth.d2App + 2 * (apparentDepth.d3App + d4AppInUse))) * x * scanVar.deltaFreqMax;
freqDeflect(2)= systemVar.diffractionMode * 2 * d4AppInUse / (2 * d4AppInUse + apparentDepth.d3App + apparentDepth.d2App) * y * scanVar.deltaFreqMax;
freqDeflect(3)= systemVar.diffractionMode * d4AppInUse / (d4AppInUse+apparentDepth.d3App) * x * scanVar.deltaFreqMax;
freqDeflect(4)= -systemVar.diffractionMode * y * scanVar.deltaFreqMax;
end