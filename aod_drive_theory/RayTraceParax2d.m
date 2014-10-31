function RayTraceParax2d()

lambda = 800e-9;
numberOfAods = 2;
V = [613, 750];
l = [10e-2, 1];
vx = 0;
linearChirps = [FirstChirp, SecondChirp];
xRange = -0.01:0.0004:0;
tRange = -1.61e-5 * (0:200);
L = 0.01;
aodEnds = [-1; 0] * L;

for t = tRange
    CalculateRayEnds(t, xRange, aodEnds);
    pause(0.1)
    clf
end

    function [xEnd,aberrationVar] = CalculateRayEnds(t, xRange, aodEnds)
        defangle = cell(1,numberOfAods); 
        x = cell(1,numberOfAods+1); 
        z = cell(1,numberOfAods+1); 
        x{1} = xRange;
        z{1} = 0.1;
        
        line([x{1};x{1}],[0;z{1}]);
        hold on
        for n = 1:numberOfAods
            x{n+1} = GetNextXyPositions(n);
            z{n+1} = z{n} + l(n);
            subplot(2,2,n)
            scatter(xRange,NthDeflectionAngle(xRange, n) * V(n) / lambda, 'r.');
            xlim([-L, 0])
            ylim([-6e7,6e7])
            xlabel('x')
            ylabel('freq')            
            subplot(2,2,[3,4])
            line([x{n};x{n+1}],[z{n};z{n+1}])
            line(aodEnds,repmat([z{n}],1,2))
            xlim([-L, 0])
            xlabel('x')
            ylabel('z')
        end
        hold off
        
        xEnd = x{end} - vx * t;
        aberrationVar = var(xEnd);
             
        function xNext = GetNextXyPositions(n)
            defangle = NthDeflectionAngle(x{1}, 1);
            for m = 2:n % does this hit if n < 2
                defangle = defangle + NthDeflectionAngle(x{m}, m);
            end
            xCurrent = x{n};
            xNext = xCurrent + defangle * l(n);
            xNext(xCurrent > aodEnds(2)) = inf;
            xNext(xCurrent < aodEnds(1)) = inf;
        end   
        
        function fnangle = NthDeflectionAngle(xin, n) 
            phase = t - xin / V(n);
            fnangle = linearChirps(n) * ShiftedPhase(phase, n) * lambda/V(n);
        end
    end

    function c = FirstChirp()
        c = V(1)^2/lambda * (V(2) - vx)/(l(1)*V(2) - l(2)*V(1) + l(2)*V(2) - l(1)*vx);
    end

    function c = SecondChirp()
        c = V(2)^2/lambda * (V(1) - vx)/(l(2)*V(1) - l(2)*V(2));
    end

    function val = ShiftedPhase(phase, n)
        adjustment = L * (1 / V(1) - 1 / V(2)); % positive
        adjustedPhase = phase;
        
        if n == 2
            adjustedPhase = adjustedPhase + adjustment;
        end
        
        factor = L / V(1);
        num = floor(adjustedPhase / factor); % add 1 to account for adjustment needed for second aod offset
        
        if n == 2
            adjustedPhase = adjustedPhase - adjustment;
        end
        
        val = adjustedPhase - num * factor;
    end

    function a = ConstAngle(phase, n)
    	factor = L / V(n);
        num = ceil(phase / factor);
        if n == 2
            a = num * L / l(n);
        else
            a = 0;
        end
    end
end
