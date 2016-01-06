function drive_equation_analysis(focalLength, correction, angleSamples, aolFunction)

lambda = 800e-9;
V = 613;
vx = 0;
vy = 0;
[t, x1, y1, k1] = GetSamplingRanges(lambda, angleSamples);

[x, y, z, k] = CalculateRayEnds(t, x1, y1, k1, focalLength, correction, aolFunction);
[x2, y2, k2] = Relay(x{end-1}, y{end-1}, k{end}, 0.8, 0.0);
[x3, y3, k3] = Objective(x2, y2, k2);
[xF, yF, d] = GetFocusedXy(k3, x3, y3);

%PlotAolFocus(x,y,z,k,t)
PlotAol(x,y,z,k)
set(gca,'FontSize',20)
PlotObjFocus(x3, y3, xF, yF, k3, d, t)
set(gca,'FontSize',20)
PlotAroundZeroPlane(k3,x3,y3,5e-4)
set(gca,'FontSize',20)
%TraceRays(x, y, z, k{end})
%TraceRays2D(x, y, z, k{end})

    function PlotObjFocus(x3, y3, xF, yF, k3, d, t)
        %TraceRays2D({x3, xF}, {y3, yF}, {0, d}, k3);
        XyScatterPsf(xF, yF, 0, 0, t);
        
        zIntercepts = x3 .* k3(3,:)./k3(1,:); % projected into xz plane
        zIntercepts = zIntercepts(~isnan(zIntercepts));
        
        zFocalPoint = mean(zIntercepts)
        zStd = std(zIntercepts - zFocalPoint)
    end

    function PlotAolFocus(x,y,z,k,t)
        TraceRays({x{end-1}, x{end}}, {y{end-1}, y{end}}, {z{end-1}, z{end}}, k{end});
        XyScatterPsf(x{end}, y{end}, 0, 0, t);
    end
        
    function PlotAol(x,y,z,k)
        TraceRays(x, y, z, k{end});
    end

    function ViewCrossSections(x,y,x2,y2,x3,y3,t)
        XyScatterPsf(x{end-1}, y{end-1}, 0, 0, t);
        XyScatterPsf(x2, y2, 0, 0, t);
        XyScatterPsf(x3, y3, 0, 0, t);
    end

    function PlotAroundZeroPlane(k3,x3,y3,D)
        [xS, yS] = GetNextXy(k3, x3, y3, -D);
        [xE, yE] = GetNextXy(k3, x3, y3, D);
        TraceRays2D({xS, xE}, {yS, yE}, {-D, D}, k3);
    end

    function [x, y, z, k] = CalculateRayEnds(t, x1, y1, k1, focalLength, cubicChirp, aolFunction)
        [aodDirectionVectors, aodSpacing, chirpFactor] = aolFunction(focalLength, 0);
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

function [xOut, yOut, d] = GetFocusedXy(k, x, y)
   
    function v = f(d) 
        [xTemp, yTemp] = GetNextXy(k, x, y, d);
        v = var(xTemp) + var(yTemp);
        v = v * 1e10;
    end
    
    d = fminsearch(@f, 0);
    [xOut, yOut] = GetNextXy(k, x, y, d);
end

function [x2, y2, k2] = Relay(x, y, k, mag, fac)
    x2 = -x*mag - fac*k(1,:)./k(3,:);
    y2 = -y*mag - fac*k(2,:)./k(3,:);
    k2 = k;
    k2(1:2,:) = -k(1:2,:)./mag;
end

function [x2, y2, k2] = Objective(x, y, k)
    r = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    phi = atan2(k(2,:), k(1,:));
    angleIn = atan(sqrt(k(1,:).^2 + k(2,:).^2) ./ k(3,:));
    
    focalLength = 0.008;
    r2 = focalLength*angleIn;
    angleOut = -asin(r/focalLength);
    
    k2 = ConvertAngleToWavevec(angleOut, theta);
    x2 = r2 .* cos(phi);
    y2 = r2 .* sin(phi);
    if ~isreal(k2)
        error('make relay mag smaller')
    end
end

function XyScatterPsf(x, y, xVel, yVel, t)
    xEff = x - xVel * t;
    yEff = y - yVel * t;
    xE = xEff(:);
    yE = yEff(:);
    figure();
    scatter(xE, yE, '.');
    axis([-1.5e-6, 1.5e-6, -1.5e-6, 1.5e-6])
    r = 2e-5;
    %axis([-r r -r r])
    xlabel('x')
    ylabel('y')
    xyStd = sqrt(var(xEff) + var(yEff)) % sum x and y vars
end

function TraceRays2D(x, y, z, kEnd)
    numberOfAods = length(x) - 1;    
    figure();
    hold on
    for n = 1:numberOfAods
        line([x{n};x{n+1}], [z{n};z{n+1}], 'color', 'red');
    end
    extra = 0.3 * (z{n+1} - z{n});
    [xNext, ~] = GetNextXy(kEnd, x{n+1}, y{n+1}, extra);
    line([x{n+1};xNext],[z{n+1};z{n+1}+extra], 'color', 'red');
    hold off
    alpha(0.1)
    axis([-5e-6, 5e-6, -15e-6, 0])
end

function TraceRays(x, y, z, kEnd)
    numberOfAods = length(x) - 1;
    figure();
    hold on
    for n = 1:numberOfAods
        line([x{n};x{n+1}],[y{n};y{n+1}],[z{n};z{n+1}]);
        fill3([max(x{n}) max(x{n}) min(x{n}) min(x{n})],...
                [max(y{n}) min(y{n}) min(y{n}) max(y{n})],...
                repmat([z{n}],1,4),repmat([z{n}],1,4))
    end
    extra = 0.3 * (z{n+1} - z{n});
    [xNext, yNext] = GetNextXy(kEnd, x{n+1}, y{n+1}, extra);
    line([x{n+1};xNext], [y{n+1};yNext], [z{n+1};z{n+1}+extra]);
    hold off
    alpha(0.1)
end

function [t,x1,y1,k1] = GetSamplingRanges(lambda, angleSamples)
    thetaRange = linspace(0,2*pi,angleSamples);
    thetaRange = thetaRange(2:length(thetaRange));
    rRange = linspace(-8,8,9) * 1e-3;
    tRange = 0;%[-1 0 1] * 1e-6;
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

