function [ eff ] = aol_analysis( rayBundle, numAodsToOptimize, scanSpeed )

effEachRay = prod(rayBundle.eff(1:numAodsToOptimize,:,:),1);
numOfDrives = rayBundle.numOfDrives;
numOfPositions = rayBundle.numOfPositions;

if sum(abs(scanSpeed)) > 0
    eff = AnalyseScanningMode(effEachRay,numOfDrives,numOfPositions);
else
    eff = AnalysePointingMode(effEachRay,numOfDrives);
end

    function [ effOverSimultaneousRays ] = AnalyseScanningMode(effEachRay,numOfDrives,numOfPositions) % want to find the efficiency by end location (time dependent) for each drive and perturbation
        effEachRay = reshape(effEachRay,size(effEachRay,2)/numOfDrives/numOfPositions,numOfPositions,numOfDrives,size(effEachRay,3));
        effOverSimultaneousRays = squeeze(mean(effEachRay,2)); % average rays for each [times x drive x perturbation]
    end

    function [ effForDrivePerturb ] = AnalysePointingMode(effEachRay,numOfDrives) % want to find the efficiency by end location for each drive and perturbation
        effEachRay = reshape(effEachRay,size(effEachRay,2)/numOfDrives,numOfDrives,size(effEachRay,3));
        effForDrivePerturb = mean(effEachRay,1); % average rays for each [drive x perturbation]
    end
end

