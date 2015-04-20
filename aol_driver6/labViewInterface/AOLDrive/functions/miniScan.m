function [tRamp, A, B, C] = miniScan( scanVar,systemVar,xyzNormalised)

% Modify the code and the functions used inside to make sure that it can
% do full frame imaging. if we have just one row of xyz start coordinates
% just do a full frame imaging and if the array has more than one row then
% do the miniScan Imaging.

%%Following Code runs for all the miniScans of interestxyzLength = xyzNormalised.imageStopNormalised - ...

xyzLength = xyzNormalised.imageStopNormalised - ...
    xyzNormalised.imageStartNormalised;

numberOfMiniScans = size(xyzNormalised.imageStartNormalised);

xyLength2 = sqrt(xyzLength(:,1).^2 + xyzLength(:,2).^2);
xyzLength2 = sqrt(xyzLength(:,1).^2 + xyzLength(:,2).^2 + xyzLength(:,3).^2);

numxyzvoxels = round(xyzLength2*scanVar.mVoxels/2) + 1;
dataCollectionTime = ceil ((numxyzvoxels * ...
    scanVar.dwellTime)/systemVar.dataClock)*systemVar.dataClock;
halfpixel = 1/scanVar.mVoxels;

costheta = xyzLength(:,1)./xyLength2;
sintheta = xyzLength(:,2)./xyLength2 ;
sinphi = xyzLength(:,3)./xyzLength2;

deltax = costheta*halfpixel;
deltay = sintheta*halfpixel;
deltaz = halfpixel*sinphi; 
deltaxyz = [deltax deltay deltaz];

[ xyzNormalised,xyzMid,xyzLength ] = ...
      CorrectxyzNormalized( xyzNormalised,deltaxyz, xyzLength);
 [apparentDepth] = depthParameters(scanVar,systemVar, xyzNormalised, xyzLength,dataCollectionTime,xyzMid);

tScan = repmat((systemVar.aodFill + scanVar.dataCollectionTime),1,4);

tStart = systemVar.aodFill-tScan./2;
tStop = tScan./2;
tMid = systemVar.aodFill/2;

[ aStart ] = aParameters( scanVar,systemVar,apparentDepth, ...
    xyzNormalised.imageStartNormalised(:,3)); % Calculates aSart using the same equations as aRamp
[ aStop ] = aParameters( scanVar,systemVar,apparentDepth, ...
    xyzNormalised.imageStopNormalised(:,3)); % Calculates aStop using the same equations as aRamp
[ aMid ] = aParameters( scanVar,systemVar,apparentDepth, ...
    xyzMid(:,3));

[ aP,bP,cP ] = ABCp( scanVar,systemVar,apparentDepth,aStart,aStop,aMid,xyzMid,tStart,tStop,tMid,numberOfMiniScans); %Calculates Ap,Bp,Cp values from aStart,aStop and apparentDepth
[ tRamp, A, B, C] = scaleTrampABC( aP,bP,cP,tScan,systemVar ); % Scales the Ap,Bp, Cp to A,B,C for FPGA
end