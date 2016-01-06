classdef WavefrontPair
    properties
        lambda
        V
        L
        F0
    end
    methods
        function s = WavefrontPair(lambda, V, L, F0)            
            s.lambda = lambda;
            s.V = V;
            s.L = L;
            s.F0 = F0;
        end
        
        function [C, Z0, v, X0] = apply_start_stop_constraints(s, x_start, x_stop, z_start, z_stop, t_half, D)
            sum_z_reciprocals = z_stop.^-1 + z_start.^-1 + (z_stop + z_start == 0) * 1e-6; % avoid inf
            Z0 = 2 ./ (sum_z_reciprocals - 16 .* s.V.^2 .* t_half.^2 .* D);
            C = Z0 ./ 2 ./ t_half .* (z_stop.^-1 - z_start.^-1);
            
            D_alt = 8 .* s.V.^2 .* Z0 .* D;
            v = ( (1 + C.*t_half + D_alt.*t_half.^2).*x_stop - (1 - C.*t_half + D_alt.*t_half.^2).*x_start ) ./ 2 ./ t_half;
            X0 =( (1 + C.*t_half + D_alt.*t_half.^2).*x_stop + (1 - C.*t_half + D_alt.*t_half.^2).*x_start ) ./ 2;
        end
        
        function [w1, w2] = calc_wavefronts_at_ref(s, x_start, x_stop, z_start, z_stop, t_half, D)
            % assumes -1 mode
            [C, Z0, v, X0] = s.apply_start_stop_constraints(x_start, x_stop, z_start, z_stop, t_half, D);
            
            if 0 && any(D)
                if any(C) || any(v)
                    w1 = WavefrontPair.quartic1(C,D,s.F0,s.L,s.V,Z0,s.lambda,v);
                    w2 = WavefrontPair.quartic2(C,D,s.F0,s.V,X0,Z0,s.lambda,v);
                else
                    w1 = WavefrontPair.quartic1_simple(D,s.F0,s.L,s.V,Z0,s.lambda);
                    w2 = WavefrontPair.quartic2_simple(D,s.F0,s.V,X0,Z0,s.lambda);
                end
            elseif any(C)
                w1 = WavefrontPair.cubic1(C,s.F0,s.L,s.V,Z0,s.lambda,v);
                w2 = WavefrontPair.cubic2(C,s.F0,s.V,X0,Z0,s.lambda,v);
            else 
                w1 = WavefrontPair.quad1(s.F0,s.V,Z0,s.lambda,v,s.L);
                w2 = WavefrontPair.quad2(s.F0,s.V,Z0,s.lambda,v,X0);
            end
        end
    end
    
    methods (Static)
        w = quad1(F0,V,Z0,lambda,v,L)
        function w = quad2(F0,V,Z0,lambda,v,X0)
            w = [-X0./Z0-(F0.*lambda)./V;...
                (1./2)./Z0-(v.*(1./2))./(V.*Z0);...
                0*X0;...
                0*X0];
        end
                
        w = cubic1(C,F0,L,V,Z0,lambda,v)
        function w = cubic2(C,F0,V,X0,Z0,lambda,v)
            w = [-X0./Z0-(F0.*lambda)./V;...
                (1.0./2.0)./Z0-(v.*(1.0./2.0))./(V.*Z0);...
                (C.*(1.0./2.0))./(V.*Z0);...
                0*X0];
        end
        
        w = quartic1(C,D,F0,L,V,Z0,lambda,v)
        function w = quartic2(C,D,F0,V,X0,Z0,lambda,v)
            w = [-X0./Z0-(F0.*lambda)./V;...
                (1.0./2.0)./Z0-(v.*(1.0./2.0))./(V.*Z0);...
                (C.*(1.0./2.0))./(V.*Z0);...
                D.*1.2e1];
        end
        
        w = quartic1_simple(D,F0,L,V,Z0,lambda)
        function w = quartic2_simple(D,F0,V,X0,Z0,lambda)
            w = [-X0./Z0-(F0.*lambda)./V;...
                (1.0./2.0)./Z0;...
                0*D;...
                D.*1.2e1];
        end
    end
end