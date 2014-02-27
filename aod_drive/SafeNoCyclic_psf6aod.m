function psf6aod()

lambda = 800e-9;
V = 613;
[t,x1,y1] = GetSamplingRanges();
xVel = 00;
yVel = 00;
focusDeviation = -0.000;

[aodDirectionVectors, aodL, chirpFactor] = Aod6pairedAS();

tic
[xy,aberrationVar] = CalculateRayEnds();
toc

figure();
hist3(xy);
xlabel('x')
ylabel('y')

    function [xy,aberrationVar] = CalculateRayEnds()
        numberOfAods = length(aodL);
        linearChirps = arrayfun(@NthAodChirp, 1:length(aodL));
                                                                                   
        x = cell(1,numberOfAods+1); 
        y = cell(1,numberOfAods+1); 
        z = cell(1,numberOfAods+1); 
        x{1} = x1;
        y{1} = y1;
        z{1} = 0;
        figure();
        hold on
        for n = 1:numberOfAods
            [x{n+1},y{n+1}] = GetNplus1xy(n);  % Stores the (n+1)th [x,y]
            z{n+1} = z{n} + aodL(n);
            line([x{n};x{n+1}],[y{n};y{n+1}],[z{n};z{n+1}]);
            fill3([max(x{n}) max(x{n}) min(x{n}) min(x{n})],[max(y{n}) min(y{n}) min(y{n}) max(y{n})],repmat([z{n}],1,4),repmat([z{n}],1,4))
        end
        hold off
        alpha(0.1)
        xEnd = x{end} - xVel * t;
        yEnd = y{end} - yVel * t;
        xy = [xEnd(:) yEnd(:)];
        aberrationVar = var(xEnd) + var(yEnd); % sum x and y vars
             
        function [xNplus1,yNplus1] = GetNplus1xy(n)
            angle = NthDeflectionAngle(1);
            for m = 2:n
                angle = angle + NthDeflectionAngle(m);
            end
            distance = aodL(n);
            if n == length(aodL)
                distance = distance + focusDeviation;
            end
            xNplus1 = x{n} + angle(1,:) * distance;
            yNplus1 = y{n} + angle(2,:) * distance;
        end   
        
        function fnangle = NthDeflectionAngle(n) 
            kn = aodDirectionVectors{n};
            kx = kn(1);
            ky = kn(2);
            phase = t - ( x{n}*kx + y{n}*ky ) / V;
            fnangle = aodDirectionVectors{n} * linearChirps(n) * phase * lambda/V;
        end
        
        function chirp = NthAodChirp(n)
            chirp = V*V/lambda * chirpFactor(n);
        end
    end

    function [aodDirectionVectors, aodL, chirpFactor] = Aod4()
        aodDirectionVectors = {[1;0], [0;1]};
        aodDirectionVectors = [aodDirectionVectors, cellfun(@(x) x*(-1), aodDirectionVectors,'UniformOutput', false)];
        aodL = [5, 5, 5, 5e1];
        l1 = aodL(1);
        l2 = aodL(2);
        l3 = aodL(3);
        l4 = aodL(4);
        chirpFactor = [ 1/(l1 + l2 + 2*l3 + 2*l4)...
                        1/(l2 + l3 + 2*l4)...
                           1/(2*(l3 + l4))...
                                  1/(2*l4)];
    end

    function [aodDirectionVectors, aodL, chirpFactor] = Aod6pairedT()
        aodDirectionVectors = {[1;0], [-1;0], [1;sqrt(3)]/2, -[1;sqrt(3)]/2, [-1;sqrt(3)]/2, -[-1;sqrt(3)]/2};
        aodL = [5, 5, 5 ,5 ,5 , 5e0];
        L = aodL(1);
        f = aodL(6);
        chirpFactor =  [(96*L^2 + 96*L*f + 24*f^2)/(422*L^3 + 612*L^2*f + 270*L*f^2 + 36*f^3)...
        -(2*L*(14*L + 6*f))/(326*L^3 + 516*L^2*f + 246*L*f^2 + 36*f^3)...
                                     -(4*L)/(30*L^2 + 42*L*f + 12*f^2)...
                               (16*L + 8*f)/(34*L^2 + 42*L*f + 12*f^2)...
                                   (10*L + 4*f)/((L + f)*(13*L + 6*f))...
                                                                     0];
    end
    function [aodDirectionVectors, aodL, chirpFactor] = Aod6pairedS()
        aodDirectionVectors = {[1;0], [-1;0], [1;sqrt(3)]/2, -[1;sqrt(3)]/2, [-1;sqrt(3)]/2, -[-1;sqrt(3)]/2};
        aodL = [5, 5, 5 ,5 ,5 , 5e0];
        L = aodL(1);
        f = aodL(6);
        chirpFactor =  [(144*L^2 + 168*L*f + 48*f^2)/((2900*L^3)/3 + 1912*L^2*f + 984*L*f^2 + 144*f^3)...
 ((56*L^2)/3 + 120*L*f + 48*f^2)/(2*((1234*L^3)/3 + 872*L^2*f + 468*L*f^2 + 72*f^3))...
                               (2*((2*L)/3 + 4*f))/(3*(14*L^2 + (80*L*f)/3 + 8*f^2))...
                                     (2*(6*L + 4*f))/((122*L^2)/3 + 72*L*f + 24*f^2)...
                                          (14*L + 4*f)/(((2*L)/3 + f)*(26*L + 12*f))...
                                                                             1/(3*f)];
    end
    function [aodDirectionVectors, aodL, chirpFactor] = Aod6pairedAT()
        aodDirectionVectors = {[1;0], [-1;0], -[1;sqrt(3)]/2, [1;sqrt(3)]/2, [-1;sqrt(3)]/2, -[-1;sqrt(3)]/2};
        aodL = [5, 5, 5 ,5 ,5 , 5e0];
        L = aodL(1);
        f = aodL(6);
        chirpFactor =  [(336*L^2 + 368*L*f + 96*f^2)/(1224*L^3 + 2128*L^2*f + 1032*L*f^2 + 144*f^3)...
          -(192*L^2 + 96*f*L)/(2*(444*L^3 + 880*L^2*f + 468*L*f^2 + 72*f^3))...
                              (2*(12*L + 8*f))/(3*(20*L^2 + 28*L*f + 8*f^2))...
                                          -(16*L)/(36*L^2 + 68*L*f + 24*f^2)...
                                        (20*L + 8*f)/((L + f)*(22*L + 12*f))...
                                                                           0];
    end
    function [aodDirectionVectors, aodL, chirpFactor] = Aod6pairedAS()
        aodDirectionVectors = {[1;0], [-1;0], -[1;sqrt(3)]/2, [1;sqrt(3)]/2, [-1;sqrt(3)]/2, -[-1;sqrt(3)]/2};
        aodL = [5, 5, 5 ,5 ,5 , 5e0];
        L = aodL(1);
        f = aodL(6);
        chirpFactor =  [((392*L^2)/3 + 168*L*f + 48*f^2)/((2308*L^3)/3 + 1704*L^2*f + 936*L*f^2 + 144*f^3)...
    (- 16*L^2 + 88*L*f + 48*f^2)/(2*((958*L^3)/3 + 768*L^2*f + 444*L*f^2 + 72*f^3))...
                             (2*((14*L)/3 + 4*f))/(3*((38*L^2)/3 + 24*L*f + 8*f^2))...
                                -(2*((2*L)/3 - 4*f))/((86*L^2)/3 + 64*L*f + 24*f^2)...
                                     ((38*L)/3 + 4*f)/(((2*L)/3 + f)*(22*L + 12*f))...
                                                                            1/(3*f)];
    end
    function [aodDirectionVectors, aodL, chirpFactor] = Aod6pairedASscan()
        aodDirectionVectors = {[1;0], [-1;0], -[1;-sqrt(3)]/2, [1;-sqrt(3)]/2, -[1;sqrt(3)]/2, [1;sqrt(3)]/2};
        aodL = [5, 5, 5 ,5 ,5 , 5e1];
        L = aodL(1);
        f = aodL(6);
        vx = 0;
        vy = 2;
        xVel = vx*V;
        yVel = vy*V;
        chirpFactor =  [(- 27*L^2*vx^2 - 10*3^(1/2)*L^2*vx*vy + 170*L^2*vx + 19*L^2*vy^2 + 130*3^(1/2)*L^2*vy + (392*L^2)/3 - 18*L*f*vx^2 + 228*L*f*vx + 6*L*f*vy^2 + 116*3^(1/2)*L*f*vy + 168*L*f + 72*f^2*vx + 24*3^(1/2)*f^2*vy + 48*f^2)/(- 27*L^3*vx^2 - 10*3^(1/2)*L^3*vx*vy + 32*L^3*vx + 19*L^3*vy^2 + (248*3^(1/2)*L^3*vy)/3 + (2308*L^3)/3 - 18*L^2*f*vx^2 + 72*L^2*f*vx + 6*L^2*f*vy^2 + 128*3^(1/2)*L^2*f*vy + 1704*L^2*f + 36*L*f^2*vx + 36*3^(1/2)*L*f^2*vy + 936*L*f^2 + 144*f^3)...
            -(- 27*L^2*vx^2 - 10*3^(1/2)*L^2*vx*vy + 206*L^2*vx + 19*L^2*vy^2 + (466*3^(1/2)*L^2*vy)/3 + 16*L^2 - 18*L*f*vx^2 + 252*L*f*vx + 6*L*f*vy^2 + 124*3^(1/2)*L*f*vy - 88*L*f + 72*f^2*vx + 24*3^(1/2)*f^2*vy - 48*f^2)/(2*(444*L*f^2 + 768*L^2*f - 69*L^3*vx + (958*L^3)/3 + 72*f^3 - (71*3^(1/2)*L^3*vy)/3 - 18*L*f^2*vx - 78*L^2*f*vx + 6*3^(1/2)*L*f^2*vy + 6*3^(1/2)*L^2*f*vy))...
            (2*((14*L)/3 + 4*f + 3*L*vx + (17*3^(1/2)*L*vy)/3 + 4*3^(1/2)*f*vy))/(3*(24*L*f - L^2*vx + (38*L^2)/3 + 8*f^2 - 2*L*f*vx + (5*3^(1/2)*L^2*vy)/3 + 2*3^(1/2)*L*f*vy))...
            (2*((2*L)/3 - 4*f + 5*L*vx + (23*3^(1/2)*L*vy)/3 + 4*3^(1/2)*f*vy))/(9*L^2*vx - 64*L*f - (86*L^2)/3 - 24*f^2 + 6*L*f*vx + (19*3^(1/2)*L^2*vy)/3 + 2*3^(1/2)*L*f*vy)...
            ((38*L)/3 + 4*f - L*vx + (5*3^(1/2)*L*vy)/3)/(((2*L)/3 + f)*(22*L + 12*f - 3*L*vx + 3^(1/2)*L*vy))...
            1/(3*f)];
    end  
end

function [t,x1,y1] = GetSamplingRanges()
    xRange = normrnd(0,10,1,10);
    yRange = xRange;
    tRange = 0;%:0.010:0.100;
    [x1temp,y1temp] = meshgrid(xRange,yRange);
    [x1temp,ttemp] = meshgrid(x1temp(:),tRange);
    [y1temp,~] = meshgrid(y1temp(:),tRange);
    x1 = x1temp(:)';
    y1 = y1temp(:)';
    t = ttemp(:)';
end

