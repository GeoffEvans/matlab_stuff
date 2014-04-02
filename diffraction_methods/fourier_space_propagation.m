function run()

rect = @(x) x <= 0.5 & x > -0.5;
delta = 200;
Y = 10;
func = @(y) rect((y - 2*delta)/Y) + rect((y - delta)/Y) + rect(y/Y) + rect((y + delta)/Y) + rect((y + 2*delta)/Y);
func = @(y) rect(y/10);

figure()

fourier_space_propagation( func, 0.17, -600, 200, 1, 100)

    function [] = fourier_space_propagation( ScalarFunction, xResolution, xStart, xEnd, k, zMax)
        
        xNumSamples = 2^ceil(log2((xEnd - xStart)/xResolution)); % round to even num
        xLength = xResolution * xNumSamples;
        xEnd = xStart + xLength;
        
        [ kx ,ft ] = dft( ScalarFunction, xResolution, xStart, xNumSamples );
        
        zVals = zMax;
        ftPropagated = PropagateFt(ft, kx, zVals);
        
        xPropagated = zeros(length(zVals),xNumSamples);
        for n = 1:length(zVals)
            xPropagated(n,:) = ifft(ftPropagated(n,:));
        end
        
        xVals = 0:xResolution:(xLength-xResolution);
        xVals = xVals - xLength * (xVals >= xEnd);
        [xGrid,zGrid] = meshgrid(xVals,zVals);
         
        [~,ind] = sortrows(xGrid');
        xGrid = xGrid(:,ind');
        xPropagated = xPropagated(:,ind');
        plot(xGrid,abs(xPropagated));
        %s = surf(xGrid,zGrid,abs(xPropagated));
        
        %set(s,'linestyle','none');
        xlabel('x')
        ylabel('z')
        zlabel('I')
        
        function propagatedSamples = PropagateFt(ftSamples, kx, zVals)
            kz = sqrt(k.^2 - kx.^2);
            [kzGrid,zGridL] = meshgrid(kz,zVals);
            phase = exp(1i .* kzGrid .* zGridL);
            propagatedSamples = repmat(ftSamples,length(zVals),1) .* phase;
        end
         
        
    end

end