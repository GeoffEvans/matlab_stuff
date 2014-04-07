function [eff, dTheta, dPhi] = rescattered_efficiency( iTheta, iPhi, acFreq, acTheta, acPower, transducerLength )
% Finds the diffracted beam efficiency for given incident angles and acoustic frequencies.

[effPrimary, dTheta, dPhi] = Efficiency( iTheta, iPhi );
effSecondary = Efficiency( dTheta, dPhi );
eff = effPrimary .* (1 - effSecondary);

    function [eff, dTheta, dPhi] = Efficiency( iTheta, iPhi )
    % Finds the diffracted beam efficiency for given incident angles and acoustic frequencies.

    L = transducerLength;
    H = 20e-3;
    opWavelenVac = aod3d.opWavelenVac;
    acPhi = pi/4 + 0*acTheta;
    
    [ dTheta, dPhi, phaseOffset ] = aod3d.match_phase(iTheta,iPhi,acFreq,acTheta);
    acVel = teo2.find_v_ac_min(acTheta, acPhi);
    soundAmp = sqrt( (2 * acPower) ./ (teo2.density * H * L .* acVel.^3) ); % (2.143) Xu & Stroud
    [ nOrd, nExt, ~,~] = teo2.find_n_op(iTheta);
    C = 0.12 / 4 .* nOrd.^2 .* nExt; % Xu & St. P66'

    indexShiftTerm = ( 2*pi .* C .* soundAmp .* L ./ opWavelenVac ).^2;
    phaseTerm = ( L/2 .* phaseOffset ).^2;
    eff = indexShiftTerm .* sinc( (indexShiftTerm + phaseTerm).^0.5 / pi ).^2;
    end
end

