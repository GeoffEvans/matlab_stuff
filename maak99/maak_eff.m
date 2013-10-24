function [ eff ] = maak_eff( acFreq, iAngleDeg )

%------------------------------ PARAMETERS ------------------------------%
L = 3.6e-3; % transducer length
opWavelenVac = 633e-9;
acAngle0 = 0 * pi / 180;
angularRes = 0.0005; % For matching wavevectors
%------------------------------------------------------------------------%

iAngle = iAngleDeg * pi / 180;
[ ~,dAngle,~,nOrd,nExt ] =...
    find_optimal_angles(iAngle, angularRes, acFreq, opWavelenVac);

indexShiftTerm = ( 2*pi .* L ./ opWavelenVac / 10e4 ).^2;
phaseTerm = ( (pi*L./opWavelenVac) .* (nExt.*cos(iAngle-acAngle0) - nOrd.*cos(dAngle-acAngle0)) ).^2;
eff = indexShiftTerm .* sinc( (indexShiftTerm + phaseTerm).^0.5 / pi ).^2;

    function optoElasticCoef = OptoElasticCoef()
        % This is given as a function of alpha (acoustic angle) and gamma in Maak 1999.
        % It is not even clear what gamma is but Xu&St (2.1) have only dependence acoustic angle.
        % Since we always have sound down the 110, we initially give this constant value.
        optoElasticCoef = 1;
    end
end

