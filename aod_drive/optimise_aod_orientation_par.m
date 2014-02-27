function [ opt ] = optimise_aod_orientation()

tic
    thetaVals = (0:0.02:0.1) * pi;
    phiVals = (-1:0.1:1) * pi;

    numberOfAods = 2;
    opt = [0.062831853071796, -2.513274122871835, 0.062831853071796, 0];
    opt = zeros(1,numberOfAods*2);
    sets = {thetaVals,phiVals,thetaVals,phiVals};
    cartProd = CartesianProduct(sets);
    combs = size(cartProd,1)
    maxEff = 0;
    for m = 1:combs
        [avX,avY,sigmaX,sigmaY,avEff,maxFracAngleError ] = PerformanceTimeAverage( cartProd(m,1:2:numberOfAods*2), cartProd(m,2:2:numberOfAods*2) );
        if avEff > maxEff
            maxEff = avEff;
            opt = cartProd(m,:);
        end
        if mod(m,100)==0
            m
        end
    end
toc

    function [ avX,avY,sigmaX,sigmaY,avEff,maxFracAngleError ] = PerformanceTimeAverage( theta, phi )

        microSecs = 0;%-2:1:2;
        [xEnd,yEnd,prodeff,fractionalAngleErrorMax] = PerformanceOverTime( microSecs, false );

        sigmaX = std(xEnd,1);
        sigmaY = std(yEnd,1);
        avX = mean(xEnd);
        avY = mean(yEnd);
        avEff = mean(prodeff);
        maxFracAngleError = max(fractionalAngleErrorMax);
    

        function [xEnd,yEnd,prodeff,fractionalAngleErrorMax] = PerformanceOverTime( microSecs, plotAods )
            fractionalAngleErrorMax = 0;
            acPower = 1; % Watts
            iPolAir = [1; 1i]/sqrt(2);
            V = 613;
            wavelengthVac = 800e-9;

            [aodAcDirectionVectors, aodL, linearChirps] = Aod2();
            numberOfAods = length(aodAcDirectionVectors);

            t = microSecs * 1e-6;
            x = zeros(numberOfAods+1,length(t));
            y = zeros(numberOfAods+1,length(t));
            z = zeros(1,numberOfAods+1);
            kVac = [0;0;1]*2*pi/wavelengthVac;
            k = repmat(kVac,1,length(t));
            eff = zeros(numberOfAods,length(t));
            
            for n=1:numberOfAods
                    [k, eff(n,:)] = DeflectAtAod(n,k);
                    x(n+1,:) = x(n,:) + aodL(n)*k(1,:)./k(3,:);
                    y(n+1,:) = y(n,:) + aodL(n)*k(2,:)./k(3,:);
            end
            
            if plotAods
                hold on;
                for n=1:numberOfAods
                    z(n+1) = z(n) + aodL(n);
                    fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))],repmat(z(n),1,4),repmat(z(n),1,4));
                end
                alpha(0.1)

                plot3(x,y,z);
                grid on;
                axis square;
                xlabel('x')
                ylabel('y')
                zlabel('z')
                hold off;
            end
            
            xEnd = x(end,:);
            yEnd = y(end,:);
            prodeff = prod(eff);
            
            function [kOut, eff] = DeflectAtAod(n, kIn)
                Kn = aodAcDirectionVectors{n};
                Kx = Kn(1);
                Ky = Kn(2);

                phase = t - ( x(n)*Kx + y(n)*Ky ) / V;
                localFreq = 35e6 + linearChirps(n) * phase;
                RotateToAlignAcousticWaveWithX = [1 -1 0; 1 1 0; 0 0 sqrt(2)]/sqrt(2)*[Kx, Ky 0; -Ky, Kx 0; 0 0 1];
                kInR = RotateToAlignAcousticWaveWithX*kIn; % change reference frame to have acoustic vector is <110> and z remains <001>
                kInRopt = RotationPhi(phi(n)-pi/4) * RotationTheta(theta(n))' * RotationPhi(phi(n)-pi/4)' * kInR; % rotate the crystal by theta 
                % about the axis at angle phi to the acoustic vector in the xy plane (no rotation about z), change to this new frame
                [ kOutRopt, eff, ~ ] = aod3d.aod_propagator_vector( kInRopt, ones(1,length(t)), repmat(iPolAir,1,length(t)), localFreq, repmat(acPower,1,length(t)) );
                kOutR = RotationPhi(phi(n)-pi/4) * RotationTheta(theta(n)) * RotationPhi(phi(n)-pi/4)' * kOutRopt; % change back to the unperturbed crystal frame
                kOut = transpose(RotateToAlignAcousticWaveWithX)*kOutR; % change back to lab frame

                defAngle = acos( dot(kOut,kIn) ./ (mag(kOut).*mag(kIn)) );
                isoAngle = wavelengthVac * localFreq / V;
                fractionalAngleError = abs(isoAngle./defAngle - 1);
                fractionalAngleErrorMax = fractionalAngleErrorMax + (fractionalAngleError > fractionalAngleErrorMax) .* (fractionalAngleError - fractionalAngleErrorMax);
            end

            function [aodDirectionVectors, aodL, chirp] = Aod4()
                aodDirectionVectors = {[0;1], [1;0], -[0;1], -[1;0]};
                %aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
                aodL = [5, 5, 5, 5e2];
                l1 = aodL(1);
                l2 = aodL(2);
                l3 = aodL(3);
                l4 = aodL(4);
                chirpFactor = [ 1/(l1 + l2 + 2*l3 + 2*l4)...
                    1/(l2 + l3 + 2*l4)...
                    1/(2*(l3 + l4))...
                    1/(2*l4)];
                chirp = V*V/wavelengthVac * chirpFactor;
            end
            
            function [aodDirectionVectors, aodL, chirp] = Aod2()
                aodDirectionVectors = {[1;0], -[1;0]};
                aodL = [5, 5e1];
                l1 = aodL(1);
                l2 = aodL(2);
                chirpFactor = [ 1/(l1 + 2*l2), 1/(2*l2) ];
                chirp = V*V/wavelengthVac * chirpFactor;
            end
            function mat = RotationPhi(phi) % from x towards y
                mat = [cos(phi) -sin(phi) 0;...
                        sin(phi) cos(phi) 0;...
                        0 0 1];
            end

            function mat = RotationTheta(theta) % from x to z
                mat = [cos(theta) 0 -sin(theta);... % rotate in xz plane
                            0       1   0;...
                        sin(theta) 0 cos(theta)];
            end
        end
    end

    function result = CartesianProduct(sets)
        c = cell(1, numel(sets));
        [c{:}] = ndgrid( sets{:} );
        result = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
    end
    function m = mag(v)
       m = dot(v,v);
       m = sqrt(m);
    end
end