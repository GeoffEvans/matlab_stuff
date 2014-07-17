function RayTraceParax2d()

lambda = 800e-9;
numberOfAods = 2;
V = [613, 700];
l = [5e-2, 1];
vx = 10;
linearChirps = [FirstChirp, SecondChirp];
[t, x1] = GetSamplingRanges();
CalculateRayEnds();

    function [xEnd,aberrationVar] = CalculateRayEnds()
        x = cell(1,numberOfAods+1); 
        z = cell(1,numberOfAods+1); 
        x{1} = x1;
        z{1} = 0;
        
        figure();
        hold on
        for n = 1:numberOfAods
            x{n+1} = GetNextXyPositions(n);
            z{n+1} = z{n} + l(n);
            line([x{n};x{n+1}],[z{n};z{n+1}]);
            line([max(x{n}) min(x{n})],repmat([z{n}],1,2))
        end
        hold off
        
        xEnd = x{end} - vx * t;
        aberrationVar = var(xEnd);
             
        function xNplus1 = GetNextXyPositions(n)
            angle = NthDeflectionAngle(1);
            for m = 2:n % does this hit if n < 2
                angle = angle + NthDeflectionAngle(m);
            end
            xNplus1 = x{n} + angle(1,:) * l(n);
        end   
        
        function fnangle = NthDeflectionAngle(n) 
            phase = t - x{n} / V(n);
            fnangle = linearChirps(n) * phase * lambda/V(n);
        end
    end

    function c = FirstChirp()
        c = V(1)^2/lambda * (V(2) - vx)/(l(1)*V(2) - l(2)*V(1) + l(2)*V(2) - l(1)*vx);
    end

    function c = SecondChirp()
        c = V(2)^2/lambda * (V(1) - vx)/(l(2)*V(1) - l(2)*V(2));
    end
end

function [t,x1] = GetSamplingRanges()
    xRange = normrnd(0,10,1,10);
    tRange = 0:0.010:0.100;
    [x1temp,ttemp] = meshgrid(xRange,tRange);
    x1 = x1temp(:)';
    t = ttemp(:)';
end