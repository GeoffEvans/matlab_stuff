function [ dThetaAir, dPhi, dIntensity, dPolarisation ] = aod_propagator( iThetaAir, iPhi, acFreq, iIntensity, iPolAir )
% Incident light of known polarisation is diffracted with calculated
% efficiency and emitted with a new polarisation.

    iThetaCrystal = RefractIntoCrystal(iThetaAir);
    polEfficiencies = GetPolEfficiency(iThetaCrystal,iPhi,iPolAir);
    
    [diffractionEfficiencies, dThetaCrystal, dPhi ] = aod3dFirst.rescattered_efficiency(iThetaCrystal, iPhi, acFreq);
    
    dIntensity = iIntensity .* diffractionEfficiencies .* polEfficiencies;
    [dThetaAir, dPolarisation] = RefractOutOfCrystal(dThetaCrystal, iPhi);
end

function polEff = GetPolEfficiency( iThetaCrystal, iPhi, incidentPolAir )
    % Find component of incident light that is in extraordinary mode.
    [ ~, ~, ~, pExt ] = teo2.find_n_op( iThetaCrystal );
    normalisedRefractedPol = teo2.GetPolVectorFromScalar(pExt, iPhi);
    polComponent = sum( normalisedRefractedPol .* conj(incidentPolAir) );
    polEff = polComponent .* conj(polComponent);
end

function [ angleCrystal ] = RefractIntoCrystal( angleAir )
    options = optimset('Algorithm','trust-region-reflective','display','off',...
        'JacobPattern',speye(length(angleAir)));
    angleCrystal = fsolve(@ZeroFunction, angleAir/2.26, options);

    function f = ZeroFunction( angleCrystal )
        [ ~, nExtAngle, ~, ~ ] = teo2.find_n_op( angleCrystal );
        f = nExtAngle .* sin(angleCrystal) - sin(angleAir);
    end
end

function [ thetaAir, polarisation ] = RefractOutOfCrystal( thetaCrystal, phi )
    [ nOrdAngle, ~, pOrd, ~ ] = teo2.find_n_op( thetaCrystal );
    thetaAir = asin( nOrdAngle .* sin(thetaCrystal) );
    polarisation = teo2.GetPolVectorFromScalar(pOrd,phi);
end

