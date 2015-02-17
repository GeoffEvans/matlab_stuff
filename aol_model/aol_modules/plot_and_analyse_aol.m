function [ eff ] = plot_and_analyse_aol( microSecs, xyInputMm, aolPerturbations, driveParams, numAodsToOptimize, transducerWidths, plotRays )

[ rayBundle ] = aol_model_rays( microSecs, xyInputMm, aolPerturbations, driveParams, transducerWidths );

if plotRays
   plot_aol_rays(rayBundle);
end

[ eff ] = aol_performance( rayBundle, numAodsToOptimize, driveParams.xyScanSpeed, driveParams.xyDeflectionMm * 1e-3 );
 
end

