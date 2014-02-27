function [eff, dTheta, dPhi] = rescattered_efficiency( iTheta, iPhi, acFreq )
% Finds the diffracted beam efficiency for given incident angles and acoustic frequencies.
[effPrimary, dTheta, dPhi] = Efficiency( iTheta, iPhi, acFreq );
effSecondary = Efficiency( dTheta, dPhi, acFreq );
eff = effPrimary .* (1 - effSecondary);
end

function [eff, dTheta, dPhi] = Efficiency( iTheta, iPhi, acFreq )
% Finds the diffracted beam efficiency for given incident angles and acoustic frequencies.

% BEGIN CONSTANTS
acAngle0 = 0 * pi / 180;
L = aod3dFirst.L;
power = 1.8;
H = 20e-3;
opWavelenVac = aod3dFirst.opWavelenVac;
% END CONSTANTS

[ ~,dTheta,dPhi,acInvWavelen,nOrd,nExt ] = aod3dFirst.match_phase(iTheta,iPhi,acFreq);
acVel = acFreq ./ acInvWavelen;
soundAmp = sqrt( (2 * power) ./ (teo2.density * H * L .* acVel.^3) ); % (2.143) Xu & Stroud
C = 0.12 / 4 .* nOrd.^2 .* nExt; % Xu & St. P66'

indexShiftTerm = 3.2;%( 2*pi .* C .* soundAmp .* L ./ opWavelenVac ).^2;
phaseTerm = ( (pi*L./opWavelenVac) .* (nExt.*cos(iTheta-acAngle0) - nOrd.*cos(dTheta-acAngle0)) ).^2;
eff = indexShiftTerm .* sinc( (indexShiftTerm + phaseTerm).^0.5 / pi ).^2;
end