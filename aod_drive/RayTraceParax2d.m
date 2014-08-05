function RayTraceParax2d()

lambda = 800e-9;
numberOfAods = 2;
V = [613, 700];
l = [5e-2, 1];
vx = 0;
linearChirps = [FirstChirp, SecondChirp];
xRange = 0:0.001:0.01;
tRange = -1.61e-5 * (0:10);
L = 0.01;
aodEnds = [0; 1] * L;

for t = tRange
    CalculateRayEnds(t, xRange, aodEnds);
    pause(0.3)
    clf
end

    function [xEnd,aberrationVar] = CalculateRayEnds(t, xRange, aodEnds)
        x = cell(1,numberOfAods+1); 
        z = cell(1,numberOfAods+1); 
        x{1} = xRange;
        z{1} = 0.1;
        
        line([x{1};x{1}],[0;z{1}]);
        hold on
        for n = 1:numberOfAods
            x{n+1} = GetNextXyPositions(n);
            z{n+1} = z{n} + l(n);
            line([x{n};x{n+1}],[z{n};z{n+1}]);
            line(aodEnds,repmat([z{n}],1,2))
        end
        hold off
        
        xEnd = x{end} - vx * t;
        aberrationVar = var(xEnd);
             
        function xNext = GetNextXyPositions(n)
            angle = NthDeflectionAngle(1);
            for m = 2:n % does this hit if n < 2
                angle = angle + NthDeflectionAngle(m);
            end
            xCurrent = x{n};
            xNext = xCurrent + angle(1,:) * l(n);
            xNext(xCurrent > aodEnds(2)) = inf;
            xNext(xCurrent < aodEnds(1)) = inf;
        end   
        
        function fnangle = NthDeflectionAngle(n) 
            phase = t - x{n} / V(n);
            fnangle = ConstAngle(phase, n) + linearChirps(n) * ShiftPhase(phase, n) * lambda/V(n);
        end
    end

    function c = FirstChirp()
        c = V(1)^2/lambda * (V(2) - vx)/(l(1)*V(2) - l(2)*V(1) + l(2)*V(2) - l(1)*vx);
    end

    function c = SecondChirp()
        c = V(2)^2/lambda * (V(1) - vx)/(l(2)*V(1) - l(2)*V(2));
    end

    function val = ShiftPhase(phase, n)
        factor = L / V(n);
        num = ceil(phase / factor);
        val = phase - num * factor;
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
