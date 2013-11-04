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
        
        % ----------------------------- Externally defined
                                        symbolic_christoffel_matrix()
        [ velocities ] =                find_v_ac_min( thetaMesh, phiMesh )
        [ velocities, eVecs ] =         find_v_ac_sorted( thetaMesh, phiMesh )
        [ velocities, eVecs ] =         find_v_ac( thetaMesh, phiMesh )
        [ eVals, eVecs ] =              find_v_ac_inner( thetaMesh, phiMesh )
                                        plot_acoustic_indicatrix()
        [ nOrd, nExt, pOrd, pExt ] =    find_n_op( thetaRange )
                                        plot_optical_indicatrix()
    end
    
end

