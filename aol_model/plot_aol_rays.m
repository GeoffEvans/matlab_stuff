function plot_aol_rays(rb)

    numOfAods = rb.numOfAods;
    zPlanesAod = cellfun(@(x) x(3,1), rb.aodCentres);
    numOfZPlanes = length(rb.xyz);
    x = zeros(numOfZPlanes,rb.numOfRaysPerPerturbation,rb.numOfPerturbations);
    y = zeros(numOfZPlanes,rb.numOfRaysPerPerturbation,rb.numOfPerturbations);
    z = zeros(numOfZPlanes,rb.numOfRaysPerPerturbation,rb.numOfPerturbations);
    for planeNum = 1:numOfZPlanes
        xyz = rb.xyz{planeNum};
        x(planeNum,:,:) = xyz(1,:,:);
        y(planeNum,:,:) = xyz(2,:,:);
        z(planeNum,:,:) = xyz(3,:,:);
    end

    figure()
    hold on;
    for m=1:numOfAods
        zValAsArray = repmat(zPlanesAod(m),1,4);
        fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
    end
    zValAsArray = repmat(rb.zFocusPredicted,1,4);
    fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
    zValAsArray = repmat(rb.zFocusModel,1,4);
    fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
    alpha(0.1)

    rayEffs = prod(rb.eff,1);
    maxEffPlot = max(rayEffs(:)) % output a reference for comparing ray colours
    normEff = rayEffs/maxEffPlot;
    for pertNum = 1:rb.numOfPerturbations
        for rayNum = 1:rb.numOfRaysPerPerturbation
            colorElem = normEff(1,rayNum,pertNum);
            p = plot3(x(:,rayNum,pertNum),y(:,rayNum,pertNum),z(:,rayNum,pertNum));
            set(p, 'Color', [colorElem,0,1-colorElem]);
        end
    end
    grid on;
    axis square;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off;
end