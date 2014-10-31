function aol_analysis()

lambda = 800e-9;
V = 613;
xVel = 0;
yVel = 0;
[t, x1, y1, k1] = GetSamplingRanges(lambda);

correction = 1e12;

[x, y, z, k] = CalculateRayEnds(t, x1, y1, k1, 1e0, correction, xVel/V, yVel/V);

%xShift = mean(x{end});
%yShift = mean(y{end});
%x = cellfun(@(q) q - xShift, x, 'UniformOutput', false);
%y = cellfun(@(q) q - yShift, y, 'UniformOutput', false);

TraceRays(4, x, y, z);
XyScatterPsf(x{end}, y{end}, xVel, yVel, t);
%[xOut, yOut, d] = RelayAndObj(x{end-1}, y{end-1}, k{end});
%XyScatterPsf(xOut, yOut, 0, 0, t);
%TraceRays(1, {x{end-1}, xOut}, {y{end-1}, yOut}, {0, d});

    function [x, y, z, k] = CalculateRayEnds(t, x1, y1, k1, focalLength, cubicChirp, vx, vy)
        [aodDirectionVectors, aodSpacing, chirpFactor] = aol_chirps.Aod6(focalLength);
        numberOfAods = length(aodSpacing);
        linearChirps = V*V/lambda * chirpFactor;
                                                                                   
        x = cell(1,numberOfAods+1); 
        y = cell(1,numberOfAods+1); 
        z = cell(1,numberOfAods+1);
        k = cell(1,numberOfAods+1);
        
        x{1} = x1;
        y{1} = y1;
        z{1} = 0;
        k{1} = k1;
        
        for n = 1:numberOfAods
            k{n+1} = GetNextK(k{n}, t, x{n}, y{n}, aodDirectionVectors{n}, linearChirps(n), cubicChirp);
            [x{n+1}, y{n+1}] = GetNextXy(k{n+1}, x{n}, y{n}, aodSpacing(n));
            z{n+1} = z{n} + aodSpacing(n);
        end
        
        function kOut = GetNextK(kIn, t, xIn, yIn, aodDir, linearChirp, cubicChirp)
            tEff = t - ( xIn*aodDir(1) + yIn*aodDir(2) ) / V;
            freq = linearChirp .* tEff + cubicChirp .* tEff.^3;
            kOut = kIn + [aodDir; 0] * 2*pi * freq / V;
        end       
    end
end

function [xOut, yOut] = GetNextXy(k, xIn, yIn, distance)
    xOut = xIn + k(1,:) ./ k(3,:) * distance;
    yOut = yIn + k(2,:) ./ k(3,:) * distance;
end

function [xOut, yOut, d] = RelayAndObj(x, y, k)
    [x2, y2, k2] = Relay(x, y, k, 0.5);
    [x3, y3, k3] = Objective(x2, y2, k2);
    [xOut, yOut, d] = GetFocusedXy(k3, x3, y3);
end

function [xOut, yOut, d] = GetFocusedXy(k, x, y)
   
    function v = f(d) 
        [xTemp, yTemp] = GetNextXy(k, x, y, d);
        v = var(xTemp) + var(yTemp);
        v = v * 1e10;
    end
    
    d = fminsearch(@f, 0);
    [xOut, yOut] = GetNextXy(k, x, y, d);
end

function [x2, y2, k2] = Relay(x, y, k, mag)
    x2 = -x*mag;
    y2 = -y*mag;
    k2 = k;
    k2(1:2,:) = -k(1:2,:)./mag;
end

function [x2, y2, k2] = Objective(x, y, k)
    r = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    angleIn = atan2(sqrt(k(1,:).^2 + k(2,:).^2), k(3,:));
    
    focalLength = 0.008;
    r2 = focalLength*angleIn;
    angleOut = -asin(r/focalLength);
    
    k2 = ConvertAngleToWavevec(angleOut, theta);
    x2 = r2 .* cos(theta);
    y2 = r2 .* sin(theta);
end

function XyScatterPsf(x, y, xVel, yVel, t)
    xEff = x - xVel * t;
    yEff = y - yVel * t;
    xE = xEff(:);
    yE = yEff(:);
    figure();
    scatter(xE, yE);
    xlabel('x')
    ylabel('y')
    display(var(xEff) + var(yEff)) % sum x and y vars
end

function TraceRays(numberOfAods, x, y, z)
    figure();
    hold on
    for n = 1:numberOfAods
        line([x{n};x{n+1}],[y{n};y{n+1}],[z{n};z{n+1}]);
        fill3([max(x{n}) max(x{n}) min(x{n}) min(x{n})],...
                [max(y{n}) min(y{n}) min(y{n}) max(y{n})],...
                repmat([z{n}],1,4),repmat([z{n}],1,4))
    end
    hold off
    alpha(0.1)
end

function [t,x1,y1,k1] = GetSamplingRanges(lambda)
    thetaRange = linspace(0,2*pi,8);
    thetaRange = thetaRange(2:length(thetaRange));
    rRange = linspace(-10,10,12) * 1e-3;
    tRange = 2 * 1e-6;
    xTemp = repmat(cos(thetaRange') * rRange, [1, 1, length(tRange)]);
    yTemp = repmat(sin(thetaRange') * rRange, [1, 1, length(tRange)]);
    tTemp = zeros([length(thetaRange), length(rRange), length(tRange)]);
    for n = 1:length(tRange)
        tTemp(:,:,n) = ones([length(thetaRange), length(rRange)]) * tRange(n);
    end

    x1 = xTemp(:)';
    y1 = yTemp(:)';
    t = tTemp(:)';
    k1 = [x1*0;x1*0;x1*0+2*pi/lambda];
end

function k = ConvertAngleToWavevec(ang, theta)
    k(1,:) = cos(theta) .* sin(ang);
    k(2,:) = sin(theta) .* sin(ang);
    k(3,:) = cos(ang);
end

