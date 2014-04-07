function [ dispVec, dWavevectorAir, dInten, dPol  ] = aod_propagator_vector( iWavevectorAir, iInten, iPolAir, acFreq, acPower, transducerLength )
acTheta = pi/2 + 0*iInten;
[ wavevector, iThetaAir, iPhi ] = get_angles_from_vector( iWavevectorAir );
[ dispVec, dThetaAir, dPhi, dInten, dPol ] = aod3d.aod_propagator( iThetaAir, iPhi, iInten, iPolAir, acFreq, acTheta, acPower, transducerLength );
[ dWavevectorAir ] = get_vector_from_angles( wavevector, dThetaAir, dPhi );

if abs(wavevector - 2*pi/aod3d.opWavelenVac) > 1e-6
    display('wavevector missmatch')
    abs(wavevector - 2*pi/aod3d.opWavelenVac)
end
end

