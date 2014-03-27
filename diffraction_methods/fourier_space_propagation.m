function run()

fourier_space_propagation( @Gaussian, 0.1, -26, 20, 1, 5)

    function [] = fourier_space_propagation( ScalarFunction, xResolution, xStart, xEnd, k, zMax)
        
        numSamples = 2^ceil(log2((xEnd - xStart)/xResolution)); % round to even num
        xLength = xResolution * numSamples;
        xEnd = xStart + xLength;
        
        [ kx ,ft ] = dft( ScalarFunction, xResolution, xStart, xEnd );
        
        zVals = 0:0.1:zMax;
        ftPropagated = PropagateFt(ft, kx, zVals);
        
        xPropagated = zeros(zMax+1,numSamples);
        for n = 1:length(zVals)
            xPropagated(n,:) = ifft(ftPropagated(n,:));
        end
        
        xVals = 0:xResolution:(xLength-xResolution);
        xVals = xVals - xLength * (xVals >= xEnd);
        [xGrid,zGrid] = meshgrid(xVals,zVals);
         
        [~,ind] = sortrows(xGrid');
        xGrid = xGrid(:,ind');
        xPropagated = xPropagated(:,ind');
        s = surf(xGrid,zGrid,abs(xPropagated));
        
        set(s,'linestyle','none');
        xlabel('x')
        ylabel('z')
        zlabel('I')
        
        function propagatedSamples = PropagateFt(ftSamples, kx, zVals)
            kz = sqrt(k.^2 - kx.^2);
            [kzGrid,zGrid] = meshgrid(kz,zVals);
            phase = exp(1i .* kzGrid .* zGrid);
            propagatedSamples = repmat(ftSamples,length(zVals),1) .* phase;
        end
         
        
    end

end