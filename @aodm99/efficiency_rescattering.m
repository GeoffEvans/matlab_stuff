function eff = efficiency_rescattering( iAngle, acFreq )
% Finds the diffracted beam efficiency for pairs of incident
% angles and acoustic frequencies.
% The pairs must be in sequence as two seperate horizontal arrays.

iAngleSize = size(iAngle);
acFreqSize = size(acFreq);
if iAngleSize(1) ~= 1 || acFreqSize(1) ~= 1 || iAngleSize(2) ~= acFreqSize(2)
    error('Inputs must be horizontal arrays of the same length.');
end

[effPrimary, dAngle] = aodm99.efficiency( iAngle, acFreq );
effSecondary = aodm99.efficiency( dAngle, acFreq );

eff = effPrimary .* (1 - effSecondary);

end

