function [] = fourier_space_propagation2d( ScalarFunc, xResolution, xStart, xEnd, yResolution, yStart, yEnd, k, zVal)

    xNumSamples = 2^ceil(log2((xEnd - xStart)/xResolution)); % round to even num
    xLength = xResolution * xNumSamples;
    xEnd = xStart + xLength;
    yNumSamples = 2^ceil(log2((yEnd - yStart)/yResolution)); % round to even num
    yLength = yResolution * yNumSamples;
    yEnd = yStart + yLength;

    [ kx,ky,ft ] = dft2d( ScalarFunc, xResolution, xStart, xNumSamples, yResolution, yStart, yNumSamples );

    ftPropagated = PropagateFt(ft, kx, ky, zVal);

    xPropagated = ifft2(ftPropagated);

    xVals = 0:xResolution:(xLength-xResolution);
    xVals = xVals - xLength * (xVals >= xEnd);
    yVals = 0:yResolution:(yLength-yResolution);
    yVals = yVals - yLength * (yVals >= yEnd);
    [xGrid,yGrid] = meshgrid(xVals,yVals);

    [~,xInd] = sortrows(xGrid');
    [~,yInd] = sortrows(yGrid);
    xGrid = xGrid(yInd,xInd);
    yGrid = yGrid(yInd,xInd);
    xPropagated = xPropagated(yInd,xInd);

    s = pcolor(xGrid,yGrid,abs(xPropagated));

    set(s,'linestyle','none');
    colorbar
    xlabel('x')
    ylabel('y')
    zlabel('I')

    function propagatedSamples = PropagateFt(ftSamples, kx,ky, zVal)
        kz = sqrt(k.^2 - kx.^2 - ky.^2);
        phase = exp(1i .* kz .* zVal);
        propagatedSamples = ftSamples .* phase;
    end
end