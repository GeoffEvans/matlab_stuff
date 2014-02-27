function [eff, dAngle] = efficiency( iAngle, acFreq )
% Finds the diffracted beam efficiency for pairs of incident
% angles and acoustic frequencies.
% The pairs must be in sequence as two seperate horizontal arrays.

iAngleSize = size(iAngle);
acFreqSize = size(acFreq);
if iAngleSize(1) ~= 1 || acFreqSize(1) ~= 1 || iAngleSize(2) ~= acFreqSize(2)
    error('Inputs must be horizontal arrays of the same length.');
end

acAngle0 = 0 * pi / 180;
L = aod2d.L;
power = 1.8;
H = 20e-3;
opWavelenVac = aod2d.opWavelenVac;

[ ~,dAngle,acInvWavelen,nOrd,nExt ] = aod2d.match_phase(iAngle, acFreq);

% (2.143) Xu & Stroud
acInvVel = acInvWavelen ./ acFreq;
soundAmp = sqrt( (2 * power .* acInvVel.^3) ./ (teo2.density * H * L) );

indexShiftTerm = ( 2*pi .* L ./ opWavelenVac / 2e4 ).^2;
phaseTerm = ( (pi*L./opWavelenVac) .* (nExt.*cos(iAngle-acAngle0) - nOrd.*cos(dAngle-acAngle0)) ).^2;
eff = indexShiftTerm .* sinc( (indexShiftTerm + phaseTerm).^0.5 / pi ).^2;

end

