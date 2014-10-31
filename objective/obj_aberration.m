function obj_aberration()   

ray_count = 10;
z = 1;
focalLength = 0.008;
V = 613;
wavln = 800e-9;

examine_shift()

function [x2, k2] = Forward(x, k)
    focalLength = 0.008;    
    x2 = focalLength*k;
    k2 = -asin(x/focalLength);
end

function [x, k] = Reverse(x2, k2)
    focalLength = 0.008;    
    k = x2/focalLength;
    x = -focalLength*sin(k2);
end

function PlotTripleGraph(x, first, second, third, titleStr)
    figure()
    subplot(2,2,1:2)
    plot(x, first)
    title(titleStr)
    subplot(2,2,3)
    plot(x, second)
    subplot(2,2,4)
    plot(x, third)
end

function examine_shift()
    shift = 80e-6;
    
    k2 = linspace(-1, 1, ray_count);
    x2 = shift * tan(k2);
    [x1, k1] = Reverse(x2, k2);
    p = polyfit(x1,k1,4);
    
    PlotRayTraces()  
    PlotAngleAndFreqShiftAgainstPosition()
    AnalyseAolCurvatureShift()
    
    function AnalyseAolCurvatureShift()
        % key equation is angle = p(4)x + p(2)x^3 = lambda/V * d(t-x/v)^3
        % equate coeffs:
        d = V^3 * p(2) / 2; % d as an angle rate rather than freq rate
        b = V * p(4) / 2; % approximately half for static focus
        
        t = 4e-6 * linspace(-1,1,100);
        curvature0 = b / V;
        curvature2 = 3 * d * t.^2 / V;
        % want t^2-curvature to be small compared to the t^0-curvature 
        PlotTripleGraph(t, curvature0 + curvature2, curvature0, curvature2, 'Time')
    end
    
    function PlotRayTraces()
        subplot(1,2,1)
        plot([0*x1 - z; 0*x1 + z], [x1 - z*tan(k1); x1 + z*tan(k1)])
        subplot(1,2,2)
        plot([0*x2 - z; 0*x2 + z], [x2 - z*tan(k2); x2 + z*tan(k2)])
        axis equal
    end
    
    function PlotAngleAndFreqShiftAgainstPosition()
        f = V/wavln * k1;
        PlotTripleGraph(x1, k1, p(2)*x1.^3 , p(4)*x1, 'Angles')
        PlotTripleGraph(x1, f/2+40e6, V/wavln * p(2)*x1.^3 /2, V/wavln * p(4)*x1 /2, 'Freqs')
    end

    
end

function obj_para()
    x_array = 0.01 * linspace(-1, 1, ray_count);
    k_array = zeros(1, ray_count) + 0.000;% - x_array / 100;
    
    [x2_array, k2_array] = Forward(x_array, k_array);
    
    subplot(2,2,1)
    plot([0*x_array - z; 0*x_array + z], [x_array - z*tan(k_array); x_array + z*tan(k_array)])
    
    subplot(2,2,2)
    plot([0*x2_array - z; 0*x2_array + z], [x2_array - z*tan(k2_array); x2_array + z*tan(k2_array)])

    k2_array_rev = linspace(-1, 1, ray_count);
    x2_array_rev = zeros(1, ray_count);
    
    [x_array_rev, k_array_rev] = Reverse(x2_array_rev, k2_array_rev);
    
    subplot(2,2,3)
    plot([0*x_array_rev - z; 0*x_array_rev + z], [x_array_rev - z*tan(k_array_rev); x_array_rev + z*tan(k_array_rev)])
    
    subplot(2,2,4)
    plot([0*x2_array_rev - z; 0*x2_array_rev + z], [x2_array_rev - z*tan(k2_array_rev); x2_array_rev + z*tan(k2_array_rev)])
end

function parax_contribution()   
    x = 0.01 * linspace(-1, 1, ray_count);
    
    R = 0.1;
    k_para = -x/R;
    k_proper = -atan(x/R);
    
    aberrations = k_para - k_proper;
    
    plot(x,aberrations)
    p = polyfit(x,aberrations,5)
end
end