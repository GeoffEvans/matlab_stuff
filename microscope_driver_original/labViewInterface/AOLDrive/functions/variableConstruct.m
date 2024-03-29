function [ scanVar,xyzNormalised,systemVar ] = variableConstruct...
    ( scanVariables,xyzImageStartNormalised,xyzImageStopNormalised,...
    systemVariables )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default Values of the sytem for testing code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scanVariables = [1.0000 800.0000 40.0000 4.3500 1.0000 200.0000 0.0800];
%xyzImageStartNormalised = [0 0 0];
%xyzImageStopNormalised = [0 0 0];
%systemVariables = [0 0 0 0 240e6 200e6 50e-9 -1 15 619 0.05 0.05 0.05 0.008 0.008 0.08 0.08 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constructing variables in struct for allowing functions readability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scanPrameters from userInterface

scanVar = struct('aolMode', scanVariables(1),...
    'opticalWaveLength', scanVariables(2),...
    'centreFrequency', scanVariables(3), ...
    'acceptanceAngle', scanVariables(4),...
    'zoomFactor', scanVariables(5),...
    'numberOfVoxels', scanVariables(6), ...
    'dwellTime',scanVariables(7), 'mVoxels',0, 'deltaFreqMax', 0);


%xyzNormalized coordinates for allImagingModes
xyzNormalised = struct('imageStartNormalised',xyzImageStartNormalised,...
    'imageStopNormalised',xyzImageStopNormalised,...
    'imageCentreNormalised', (xyzImageStopNormalised+xyzImageStartNormalised)/2);
  
%systemVariables, these are fixed for a particular system
systemVar = struct('skewZX',systemVariables(1),...
    'skewZY',systemVariables(2),...
    'magZX', systemVariables(3),...
    'magZY',systemVariables(4),...
    'controlClockFreq',systemVariables(5),...
    'daqClockFreq',systemVariables(6),...
    'dataTimeInterval',systemVariables(7),... 
    'diffractionMode',systemVariables(8),...
    'crystalLength',systemVariables(9),...
    'acousticVelocity',systemVariables(10),...
    'seperationOfAod1Aod2',systemVariables(11),...
    'seperationOfAod2Aod3',systemVariables(12),...
    'seperationOfAod3Aod4',systemVariables(13),...
    'thicknessOfAod1',systemVariables(14),...
    'thicknessOfAod2',systemVariables(15),...
    'thicknessOfAod3',systemVariables(16),...
    'thicknessOfAod4',systemVariables(17),...
    'pairDeflectionRatio',systemVariables(18),...
    'dataClock', 0, 'aodFill',0);


dataClockSys = (systemVar.controlClockFreq*systemVar.dataTimeInterval)/systemVar.controlClockFreq;
dataClockDaq = (systemVar.daqClockFreq*systemVar.dataTimeInterval)/systemVar.daqClockFreq;

if (dataClockSys == dataClockDaq)
    systemVar.dataClock = 2*dataClockSys;
else
   error('sysclock is not an absolute multiple of 20MHz')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute support variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xLength = xyzNormalised.imageStopNormalised(:,1) - ...
%     xyzNormalised.imageStartNormalised(:,1); % These are 1D arrays
% yLength = xyzNormalised.imageStopNormalised(:,2) - ...
%     xyzNormalised.imageStartNormalised(:,2); % These are 1D arrays
% zLength = xyzNormalised.imageStopNormalised(:,3) - ...
%     xyzNormalised.imageStartNormalised(:,3); % These are 1D arrays
% 
% xyLength = sqrt(xLength.^2 + yLength.^2); % These are 1D arrays
% xyzLength = sqrt(xLength.^2 + yLength.^2 + zLength.^2);
% 
% scanVar.mVoxels = scanVar.numberOfVoxels*scanVar.zoomFactor;
% scanVar.deltaFreqMax = scanVar.acceptanceAngle*systemVar.acousticVelocity/scanVar.opticalWaveLength;  
% systemVar.aodFill = systemVar.dataClock * floor(systemVar.crystalLength/...
%     systemVar.acousticVelocity/systemVar.dataClock);
% numxyzvoxels=round(xyzLength.*scanVar.mVoxels/2) + 1;
% scanVar.dataCollectionTime = ceil((numxyzvoxels.*scanVar.dwellTime)/systemVar.dataClock)*systemVar.dataClock;
end

