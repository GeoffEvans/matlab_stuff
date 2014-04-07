function [ displacementVector, dThetaAir, dPhi, dInten, dPol ] = aod_propagator( iThetaAir, iPhi, iInten, iPolAir, acFreq, acTheta, acPower, transducerLength )
% Incident light of known polarisation is diffracted with calculated
% efficiency and emitted with a new polarisation.

    iThetaCrystal = RefractIntoCrystal(iThetaAir);
    polEfficiencies = GetPolEfficiency(iThetaCrystal,iPhi,iPolAir);
    [diffractionEfficiencies, dThetaCrystal, dPhi ] = aod3d.rescattered_efficiency( iThetaCrystal, iPhi, acFreq, acTheta, acPower, transducerLength );
    displacementVector = GetDisplacementInCrystal(dThetaCrystal,dPhi, transducerLength);
    dInten = iInten .* diffractionEfficiencies .* polEfficiencies;
    [dThetaAir, dPol] = RefractOutOfCrystal(dThetaCrystal, iPhi);
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

function disp = GetDisplacementInCrystal(dThetaCrystal,dPhi,transducerLength)
    directionInCrystal = get_vector_from_angles(1,dThetaCrystal,dPhi);
    disp =  directionInCrystal ./ repmat(directionInCrystal(3,:),3,1) * transducerLength;
end
