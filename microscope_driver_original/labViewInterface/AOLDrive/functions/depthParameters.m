function [apparentDepth] = depthParameters(scanVar,systemVar, xyzNormalised, xyzLength,dataCollectionTime,xyzMid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates apparent depths and dphibydt
% revisit the code to use the same funciton everywhere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d1App = (systemVar.seperationOfAod1Aod2 - systemVar.thicknessOfAod1 + ...
    systemVar.thicknessOfAod1 / 2.3); 
d2App = (systemVar.seperationOfAod2Aod3 - systemVar.thicknessOfAod2 + ...
    systemVar.thicknessOfAod2 / 2.3); 
d3App = (systemVar.seperationOfAod3Aod4 - systemVar.thicknessOfAod3 + ...
    systemVar.thicknessOfAod3 / 2.3);

dtheta4bydt = (4 * scanVar.acceptanceAngle .* (xyzLength(:,1)/2) ./...
    dataCollectionTime);


d4AppMin = systemVar.acousticVelocity^2 * ((systemVar.crystalLength / ...
    systemVar.acousticVelocity)...
    + scanVar.dwellTime)/(4 * scanVar.opticalWaveLength * scanVar.deltaFreqMax);

d4App = d4AppMin/(xyzNormalised.imageCentreNormalised(3)+.0001); % get zcentre from outside.

dphibydt = (4 * scanVar.acceptanceAngle * (xyzLength(:,2) ./ 2)) ./ ...
    dataCollectionTime;
if size(xyzMid(:,3)) == 1
    d4AppMid = [d4AppMin./xyzMid(:,3)];
else
    d4AppMid = [d4AppMin./xyzMid(:,3)];% this need to be generated in the depthParameters
end

dtheta3bydt = (dtheta4bydt.*(d3App+d4AppMid))./d4AppMid;
apparentDepth = struct('d1App',d1App,'d2App',d2App,'d3App',d3App,...
    'dtheta3bydt',dtheta3bydt,'dtheta4bydt',dtheta4bydt,'d4App',d4App,...
    'd4AppMin',d4AppMin,'d4AppMid',d4AppMid,'dphibydt',dphibydt);
end

