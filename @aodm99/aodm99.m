classdef aodm99
    
    properties (Constant)
        L = 3.6e-3; % transducer length
        opWavelenVac = 633e-9;
    end
    
    methods (Static)        % Externally defined
        [ acAngleOpt,dAngleOpt, acInverseWavelenOpt, nOrdOpt, nExt ] =  match_phase( iAngle, acFreq, displayPlot )
        [ acAngleOpt,dAngleOpt, acWavelenOpt, nOrdOpt, nExt ] =  graphical_match_phase( iAngle, acFreq )
        [ eff, dAngle ] = efficiency( acFreq, iAngleDeg )
        eff = efficiency_rescattering( acFreq, iAngleDeg )
        plot_efficiency_surface()
    end
    
end

