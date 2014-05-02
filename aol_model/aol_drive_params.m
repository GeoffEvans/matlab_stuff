classdef aol_drive_params
    % Simple data structure for holding drive parameters used to calculate
    % drive frequencies
    
    properties        
        optimalBaseFreq
        xyDeflectionMm
        pairDeflectionRatio
        scanSpeed
        focalLength
    end
    
    methods
        function obj = aol_drive_params(focalLength, optimalBaseFreq, xyDeflectionMm, pairDeflectionRatio, scanSpeed)
            if numel(focalLength) * numel(optimalBaseFreq) ~= 1
                error('only allow scalar focalLength and optimalBaseFreq')
            end
            if sum(size(scanSpeed) ~= size(pairDeflectionRatio)) &&  size(xyDeflectionMm,2) ~= size(pairDeflectionRatio,2)
                error('deflections, ratios and speeds must have same widths')
            end
            if size(scanSpeed,1) ~= 1 &&  size(xyDeflectionMm,1) ~= 2
                error('inputs should be in rows')
            end            
            obj.optimalBaseFreq = optimalBaseFreq;
            obj.xyDeflectionMm = xyDeflectionMm;
            obj.pairDeflectionRatio = pairDeflectionRatio;
            obj.scanSpeed = scanSpeed;
            obj.focalLength = focalLength;
        end
    end
    
end

