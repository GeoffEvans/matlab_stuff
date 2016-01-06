classdef teo2
    
    properties (Constant)
        activityVector = 2.65e-05; %See Warner White Bonner, 87deg/mm
        density = 5990; %kg/m3
    end
    
    methods (Static)
        function f = ElasticStiffness()
            f = [[5.57 5.12 2.18; 5.12 5.57 2.18; 2.18 2.18 10.58],...
                zeros(3); zeros(3) diag([2.65 2.65 6.59])] * 1e10;
        end
        function f = RelativePermittivity()
            f = diag([2.2597^2,2.2597^2,2.4119^2]);
        end
        function f = RelativeImpermeability()
            f = inv(teo2.RelativePermittivity());
        end
        function polVector = GetPolVectorFromScalar(polScalar,phi)
            normalisation = (1 + abs(polScalar).^2).^-0.5;
            polVector = RotatePol(phi, [normalisation; polScalar .* normalisation]);
            function rotatedVectors = RotatePol( phi, matrixOfVectors )
                % polScalar assumes phi = pi/2 so need to convert to given phi.
                rotatedVectors = 0 * matrixOfVectors;
                rotatedVectors(1,:) = sin(phi) .* matrixOfVectors(1,:) + cos(phi) .* matrixOfVectors(2,:);
                rotatedVectors(2,:) = -cos(phi) .* matrixOfVectors(1,:) + sin(phi) .* matrixOfVectors(2,:);
            end
        end
        
        % ----------------------------- Externally defined
                                        symbolic_christoffel_matrix()
        [ velocities ] =                find_v_ac_min( thetaMesh, phiMesh )
        [ velocities, eVecs ] =         find_v_ac_sorted( thetaMesh, phiMesh )
        [ velocities, eVecs ] =         find_v_ac( thetaMesh, phiMesh )       
                                        plot_acoustic_indicatrix()
        [ nOrd, nExt, pOrd, pExt ] =    find_n_op( thetaRange )
                                        plot_optical_indicatrix()       
    end
    methods (Access=private, Static)
        [ eVals, eVecs ] =              find_v_ac_inner( thetaMesh, phiMesh )
    end
end

