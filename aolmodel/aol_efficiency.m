function [ eff ] = aol_efficiency( microSecs, xyInputMm, thetaPhiAodPerturbations, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed, numAodsToOptimize, plotRays )

[ rayBundle ] = aol_model( microSecs, xyInputMm, thetaPhiAodPerturbations, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed );

if plotRays
   plot_aol_rays(rayBundle);
end

[ eff ] = aol_analysis( rayBundle, numAodsToOptimize, scanSpeed );
 
end

