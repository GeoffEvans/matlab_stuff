function run()

fourier_space_propagation( @Gaussian, 3, -40, 40, 1, 5)

    function [] = fourier_space_propagation( ScalarFunction, xResolution, xStart, xEnd, k, zMax)
        
        xSize = xEnd - xStart;
        numSamples = 2^ceil(log2(xSize/xResolution));
        kResolution = 2*pi/numSamples / xResolution;
        
        xVals = xResolution * (0:numSamples-1);
        kVals = kResolution * (0:numSamples-1);
        zVals = 0:zMax;
        
        kSamples = fft(ScalarFunction(xVals));
        propagatedSamples = PropagateFt(kSamples, kVals);
        for n = 1:zMax+1
            xSamples(n,:) = ifft(propagatedSamples(n,:));
        end
        
        [x,z] = meshgrid(xVals,zVals);
        halfNs = numSamples/2;
        s = surf([x(:,halfNs:numSamples-1)-xSize,x(:,1:halfNs)],z,abs([xSamples(:,halfNs:numSamples-1),xSamples(:,1:halfNs)]));
        set(s,'linestyle','none');
        xlabel('x')
        ylabel('z')
        zlabel('I')
        
        function propagatedSamples = PropagateFt(ftSamples, kx)
            kz = sqrt(k.^2 - kx.^2);
            [kzGrid,zGrid] = meshgrid(kz,zVals);
            phase = exp(1i .* kzGrid .* zGrid);
            propagatedSamples = repmat(ftSamples,zMax+1,1) .* phase;
        end
        
    end

end