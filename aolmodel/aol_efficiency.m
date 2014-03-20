function [ eff ] = aol_efficiency( microSecs, xyInputMm, thetaPhiAodPerturbations, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed, numAodsToOptimize, plotRays )

[ rayBundle ] = aol_model( microSecs, xyInputMm, thetaPhiAodPerturbations, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, scanSpeed );

if plotRays
    PlotRays(plotRays,rayBundle);
end

[ eff ] = aol_analysis( rayBundle, numAodsToOptimize );

    function PlotRays(rayBundle)
        
        numOfAods = rayBundle.numOfAods;
        zPlanesAod = rayBundle.aodCentres(3,:);
        x = zeros(rayBundle.numOfRaysPerPerturbation,rayBundle.numOfPerturbations);
        y = zeros(rayBundle.numOfRaysPerPerturbation,rayBundle.numOfPerturbations);
        z = zeros(rayBundle.numOfRaysPerPerturbation,rayBundle.numOfPerturbations);
        for p = 1:length(rayBundle.xyz)
            xyz = rayBundle.xyz{p};
            x(p,:,:) = xyz(1,:,:);
            y(p,:,:) = xyz(2,:,:);
            z(p,:,:) = xyz(3,:,:);
        end
        
        figure()
        hold on;
        for m=1:numOfAods
            zValAsArray = repmat(zPlanesAod(m),1,4);
            fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
        end
        zValAsArray = repmat(rayBundle.zFocusPredicted,1,4);
        fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
        zValAsArray = repmat(rayBundle.zFocusModel,1,4);
        fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
        alpha(0.1)
        for q = 1:numOfPerturbations
            plot3(x(:,:,q),y(:,:,q),z(:,:,q));
        end
        set(gcf,'DefaultAxesColorOrder',cumprod(rayBundle.eff,1));
        grid on;
        axis square;
        xlabel('x')
        ylabel('y')
        zlabel('z')
        hold off;
    end
end

