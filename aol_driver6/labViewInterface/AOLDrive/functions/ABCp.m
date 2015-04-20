function [ aP,bP,cP ] = ABCp( scanVar,systemVar,apparentDepth,aStart,...
    aStop,aMid,xyzMid,tStart,tStop,tMid,numberOfMiniScans )
% Check whether this function can be used in other modes.

cP = (aStop - aStart) ./ (2 .* (tStop - tStart));
bMid = aMid - 2*cP .* tMid;
bP = bMid;
d4AppInUse = apparentDepth.d4AppMid;

for i = 1:numberOfMiniScans(1)
  d4AppInUse = apparentDepth.d4AppMid(i);  
[ freqDeflect ] = freqDeflectEqu( scanVar,systemVar,apparentDepth,d4AppInUse,xyzMid(i,1),xyzMid(i,2))

aP(i,1) = scanVar.centreFrequency + freqDeflect(1) - bP(i,1) .* tMid - cP(i,1) .* tMid.^2;
aP(i,2) = scanVar.centreFrequency + freqDeflect(2) - bP(i,2) .* tMid - cP(i,2) .* tMid.^2;
aP(i,3) = scanVar.centreFrequency + freqDeflect(3) - bP(i,3) .* tMid - cP(i,3) .* tMid.^2;
aP(i,4) = scanVar.centreFrequency + freqDeflect(4) - bP(i,4) .* tMid - cP(i,4) .* tMid.^2;

end

end

