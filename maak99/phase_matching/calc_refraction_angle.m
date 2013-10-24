function [ iAngleCrystal, nExt ] = calc_refraction_angle( iAngleAir, angularRes )

% This doesn't work - need to worry about walk off.
% We know the gradient of n increases with theta.
% Refracted angle must be smaller than "usual".
maxAngle = 1.000277 * iAngleAir;
angleRange = 0:angularRes:maxAngle;

[ ~, nExtRange, ~, ~, ~ ] = find_n_for_angle( abs(angleRange) );
d = angularRes/10;
[ ~, nExtRangeD, ~, ~, ~ ] = find_n_for_angle( abs(angleRange + d) );
dn = (nExtRangeD - nExtRange) / d;

zero = nExtRange .* angleRange + dn - 1.000277 .* iAngleAir;
minimiser = (abs(zero) == min(abs(zero)));
iAngleCrystal = angleRange(minimiser);
nExt = nExtRange(minimiser);

end

