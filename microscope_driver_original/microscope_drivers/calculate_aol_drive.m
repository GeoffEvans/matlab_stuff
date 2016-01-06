function [aodDirectionVectors, aodCentres, zFocusPredicted, aolDrives] = calculate_aol_drive(numOfAods, driveParams, aodMode)
% Returns (cell()) aodDirectionVectors; (cell[numOfAods](3,numOfDrives)) aodCentres; (numOfAods,numOfDrives) chirp; numOfAods,numOfDrives baseFreq}

optimalBaseFreq = driveParams.optimalBaseFreq;
xyDeflectionMm = driveParams.xyDeflectionMm;
pairDeflectionRatio = driveParams.pairDeflectionRatio;
xyScanSpeed = driveParams.xyScanSpeed;
focalLengths = driveParams.focalLength;
lambda = driveParams.opWavelenVac;

widthOfAod = 8e-3;
correctionDistance = widthOfAod * ( 1 - 1/2.26 );
V = teo2.find_v_ac_min(pi/2,pi/4);

[aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = Aod4(focalLengths, xyDeflectionMm, pairDeflectionRatio, optimalBaseFreq, xyScanSpeed);

chirp = V*V/lambda * chirpFactor;
[aodCentres,zFocusPredicted] = CalculateAodCentres(numOfAods,baseFreq,aodSpacing,aodDirectionVectors); % use optimalBaseFreq because this needs to be fixed over pointing fov.
aolDrives = aol_drive_freqs(baseFreq,chirp);

    function [aodCentresCell,zFocusPredicted] = CalculateAodCentres(numOfAods,baseFreq,aodSpacing,aodDirectionVectors)
        numOfDrives = size(baseFreq,2);
        
        zPlanesAod = 1 + [0, cumsum(aodSpacing)];
        zFocusPredicted = zPlanesAod(end);
        
        aodCentresCell = cell(1,numOfAods); % calc centre of AODs and the focal point using small angle approx.
        aodCentresCell{1} = [zeros(2,numOfDrives); ones(1,numOfDrives) * zPlanesAod(1)]; % first aod centre defines xy-origin
        
        for aodNumber = 2:numOfAods % start from 2 because offset on first AOD is zero
            xyCentres = XyDeflectionThroughPriorAods(aodNumber,numOfDrives,baseFreq,aodSpacing,aodDirectionVectors);
            zCentre = zPlanesAod(aodNumber);
            aodCentresCell{aodNumber} = [xyCentres; ones(1,numOfDrives) * zCentre]; 
        end
        
        function xyDeflectionCumulator = XyDeflectionThroughPriorAods(aodNumber,numOfDrives,baseFreq,aodSpacing,aodDirectionVectors)
            xyDeflectionCumulator = zeros(2,numOfDrives);
            for mthAod = 1:aodNumber-1
                effctvSpacesFromMthToAodNum = aodSpacing(mthAod:(aodNumber-1)) - correctionDistance;
                offsetDueToMthAod = sum(effctvSpacesFromMthToAodNum) .* aodDirectionVectors{mthAod} * baseFreq(mthAod,:); % note use of outer product here
                xyDeflectionCumulator = xyDeflectionCumulator + offsetDueToMthAod;
            end
            xyDeflectionCumulator = xyDeflectionCumulator * aodMode*lambda/V;
        end
    end

    function [aodDirectionVectors, aodSpacing, chirpFactor, baseFreq] = Aod4(focalLengths, xyDeflectionMm, pair_deflection_ratio, optimalBaseFreq, xyScanSpeed)
        aodDirectionVectors = {[1;0], [0;1], -[1;0], -[0;1]};
        aodSpacing = [5e-2, 5e-2, 5e-2];
        
        spacing = aodSpacing - correctionDistance;
        focus_z = focalLengths - correctionDistance; % seperate to handle many focal lengths
        
        A = xyScanSpeed / V;
        chirpFactor = [(1 + A(1,:))./(spacing(1) + spacing(2) + 2*spacing(3) + 2*focus_z + A(1,:)*spacing(1) + A(1,:)*spacing(2));...
            (1 + A(2,:))./(spacing(2) + spacing(3) + 2*focus_z + A(2,:)*spacing(2) + A(2,:)*spacing(3));... 
            (1 - A(1,:))./(2 * (spacing(3) + focus_z) );...
            (1 - A(2,:))./(2 * focus_z)] / aodMode;
        
        % for constant components, choose r to represent the ratio of ANGULAR deflection on the first of the pair to the second of the pair
        % traditionally, this means r = 1, while all on the second would be r = 0

        multiplier = V/lambda/aodMode;       
        xy_deflection = xyDeflectionMm * 1e-3;
        
        if aodMode ~= -1
            error('PDR code assumes aod mode is -1. Need to rewrite.')
        end
        
        z_x = spacing(3) + focus_z;
        L_x = sum(spacing(1:2));
        z_y = focus_z;
        L_y = sum(spacing(2:3));
        
        x0 = optimalBaseFreq(1) .* (L_x + z_x) - optimalBaseFreq(3) .* z_x;
        y0 = optimalBaseFreq(2) .* (L_y + z_y) - optimalBaseFreq(4) .* z_y;
        
        if pair_deflection_ratio == -1
            %min_pdr_x = (( multiplier .* xy_deflection(1,:) - optimalBaseFreq(1) .* (sum(spacing) + focus_z) + optimalBaseFreq(3) .* (spacing(3) + focus_z) ) ./ (optimalBaseFreq(1) - 30e6) - (spacing(3) + focus_z) ) ./ (sum(spacing) + focus_z);
            %min_pdr_y = (( multiplier .* xy_deflection(2,:) - optimalBaseFreq(2) .* (sum(spacing(2:end)) + focus_z) + optimalBaseFreq(4) .* focus_z ) ./ (optimalBaseFreq(2) - 30e6) - focus_z) ./ (sum(spacing(2:end)) + focus_z);

            % for linear LFL
            min_pdr_x = 10 + 0*focalLengths;           
            x_on_z = ( V ./ lambda .* xy_deflection(1,:) + x0 ) ./ z_x;
            idx_x = x_on_z > -18e6 .* (1 + L_x./z_x);
            min_pdr_x(idx_x) = - (18e6 + 2.15 * x_on_z(idx_x)) ./ (x_on_z(idx_x) + 18e6 * (L_x./z_x(idx_x) + 1));
            
            min_pdr_y = 10 + 0*focalLengths;
            y_on_z = ( V ./ lambda .* xy_deflection(2,:) + y0 ) ./ z_y;
            idx_y = y_on_z > -18e6 .* (1 + L_y./z_y);
            min_pdr_y(idx_y) = - (18e6 + 2.15 * y_on_z(idx_y)) ./ (y_on_z(idx_y) + 18e6 * (L_y./z_y(idx_y) + 1));

            pair_deflection_ratio_x = max(0, min_pdr_x);
            pair_deflection_ratio_y = max(0, min_pdr_y);
        else
            pair_deflection_ratio_x = pair_deflection_ratio;
            pair_deflection_ratio_y = pair_deflection_ratio;
        end
        
        dfx = ( multiplier .* xy_deflection(1,:) - x0 ) ./ ( pair_deflection_ratio_x .* (L_x + z_x) + z_x );
        dfy = ( multiplier .* xy_deflection(2,:) - y0 ) ./ ( pair_deflection_ratio_y .* (L_y + z_y) + z_y );
        
        dfx(focalLengths == inf) = 0; % for calculating coordinate reference points
        dfy(focalLengths == inf) = 0;
        
        baseFreq = [ 
            optimalBaseFreq(1) + pair_deflection_ratio_x.*dfx; 
            optimalBaseFreq(2) + pair_deflection_ratio_y.*dfy; 
            optimalBaseFreq(3) - dfx; 
            optimalBaseFreq(4) - dfy 
            ];
    end
end