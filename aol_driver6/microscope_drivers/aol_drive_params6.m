classdef aol_drive_params6
    % Simple data structure for holding drive parameters used to calculate
    % drive frequencies
    
    properties        
        optimalBaseFreq
        xyDeflectionMm
        xyScanSpeed
        focalLength
        opWavelenVac
    end
    
    methods
        function obj = aol_drive_params6(focalLength, optimalBaseFreq, xyDeflectionMm, xyScanSpeed, opWavelenVac)
            if numel(optimalBaseFreq) == 1 
                optimalBaseFreq = optimalBaseFreq * [1,1,1,1];
            end
            if numel(optimalBaseFreq) ~= 6
                error('check size of optimalBaseFreq')
            end
            if numel(opWavelenVac) ~= 1
                error('only allow scalar wavelength')
            end
            if size(xyScanSpeed,1) ~= 2 ||  size(xyDeflectionMm,1) ~= 2
                error('inputs should be in rows')
            end
            if numel(focalLength) == 1
                focalLength = repmat(focalLength, 1, size(xyScanSpeed,2));
            end
            if numel(focalLength) ~= size(xyScanSpeed,2)
                error('check focalLength size')
            end
            
            obj.optimalBaseFreq = optimalBaseFreq;
            obj.xyDeflectionMm = xyDeflectionMm;
            obj.xyScanSpeed = xyScanSpeed;
            obj.focalLength = focalLength;
            obj.opWavelenVac = opWavelenVac;
        end
    end
       
    methods(Static)
        function driveParams = MakeDriveParams(xyDef,speed,optimalBaseFreq,focalLength,opWavelenVac)
            % returns arrays of horizontal form pure stretch on speed, pure repeat on def
            [xDefStretch,xSpeedStretched] = meshgrid(xyDef(1,:),speed(1,:));
            [yDefStretch,ySpeedStretched] = meshgrid(xyDef(2,:),speed(2,:));
            xyDefStretched = [xDefStretch(:)'; yDefStretch(:)']; 
            speedStretched = [xSpeedStretched(:)'; ySpeedStretched(:)']; 
            driveParams = aol_drive_params6(focalLength, optimalBaseFreq, xyDefStretched, speedStretched, opWavelenVac);
        end
    end
    
end

