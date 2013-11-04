function plot_efficiency_surface_2()

acFreqMax = 90e6;
acFreqMin = 10e6;
range = 1:40;
acFreqRange = range/max(range)*(acFreqMax-acFreqMin) + acFreqMin;
iAngleAirDegRange = 1.6:0.05:3.2;

iAngleAirRange = iAngleAirDegRange / 180 * pi;
iAngleCrystalRange = Refract(iAngleAirRange);

[iAngleCrystalMesh, acFreqMesh] = meshgrid(iAngleCrystalRange, acFreqRange);
[iAngleAirDegMesh, ~] = meshgrid(iAngleAirDegRange, acFreqRange);
tic
efficienciesArray = aodm99.efficiency_rescattering(iAngleCrystalMesh(:)', acFreqMesh(:)');
toc

polEfficiencyFactors = GetPolEfficiency(iAngleCrystalRange,1i);
polEffMesh = repmat(polEfficiencyFactors,length(acFreqRange),1);

efficiencies = reshape(efficienciesArray, size(iAngleCrystalMesh)) .* polEffMesh;
figure();
surf(iAngleAirDegMesh, acFreqMesh / 1e6, efficiencies, 'linestyle', 'none');
xlabel('incidence angle air / degrees')
ylabel('acoustic freq / MHz')
zlabel('efficiency')
grid on;
grid minor;
axis square;

figure();
plot(iAngleAirRange, polEfficiencyFactors);
xlabel('incidence angle air / degrees')
ylabel('polarisation efficiency / %')
zlabel('efficiency')
grid on;
grid minor;
axis square;

end

function polEff = GetPolEfficiency( iAngleCrystalRange, airPol )
    [ ~, ~, ~, pExt ] = teo2.find_n_op( abs(iAngleCrystalRange) );        
    dotProd = ( 1 + pExt.*conj(airPol) );
    normalisation = sqrt( ( 1 + airPol.*conj(airPol) ) .* ( 1 + pExt.*conj(pExt) ) );
    polComponent = dotProd ./ normalisation;
        
    polEff = polComponent .* conj(polComponent);
end

function [ angleCrystal ] = Refract( angleAir )
angleCrystal = fsolve(@ZeroFunction, angleAir/2.26);

    function f = ZeroFunction( angleCrystal )
        [ ~, nExtAngle, ~, ~ ] = teo2.find_n_op( abs(angleCrystal) );
        f = nExtAngle .* sin(angleCrystal) - sin(angleAir);
    end
end