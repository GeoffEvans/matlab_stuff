function [ eff ] = simulate_aol( microSecs, xyInputMm, aolPerturbations, driveParams, numAodsToOptimize, plotRays )

[ rayBundle ] = aol_model_rays( microSecs, xyInputMm, aolPerturbations, driveParams );

if plotRays
   plot_aol_rays(rayBundle);
end

[ eff ] = aol_analysis( rayBundle, numAodsToOptimize, driveParams.scanSpeed );
 
end

