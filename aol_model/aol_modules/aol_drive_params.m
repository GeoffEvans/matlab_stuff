classdef aol_drive_params
    % Simple data structure for holding drive parameters used to calculate
    % drive frequencies
    
    properties        
        optimalBaseFreq
        xyDeflectionMm
        pairDeflectionRatio
        xyScanSpeed
        focalLength
        opWavelenVac
    end
    
    methods
        function obj = aol_drive_params(focalLength, optimalBaseFreq, xyDeflectionMm, pairDeflectionRatio, xyScanSpeed, opWavelenVac)
            if numel(optimalBaseFreq) == 1 
                optimalBaseFreq = optimalBaseFreq * [1,1,1,1];
            end
            if numel(optimalBaseFreq) ~= 4
                error('check size of optimalBaseFreq')
            end
            if numel(opWavelenVac) ~= 1
                error('only allow scalar wavelength')
            end
            if sum((size(xyScanSpeed,2) ~= size(pairDeflectionRatio,2)) + (size(xyDeflectionMm,2) ~= size(pairDeflectionRatio,2)))
                error('deflections, ratios and speeds must have same widths')
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
            obj.pairDeflectionRatio = pairDeflectionRatio;
            obj.xyScanSpeed = xyScanSpeed;
            obj.focalLength = focalLength;
            obj.opWavelenVac = opWavelenVac;
        end
    end
       
    methods(Static)
        function driveParams = MakeDriveParams(xyDef,ratio,speed,optimalBaseFreq,focalLength,opWavelenVac)
            % returns arrays of horizontal form pure stretch on speed, pure repeat on ratio
            [xDefStretch,ratioStretched,xSpeedStretched] = meshgrid(xyDef(1,:),ratio,speed(1,:));
            [yDefStretch,~,ySpeedStretched] = meshgrid(xyDef(2,:),ratio,speed(2,:));
            xyDefStretched = [xDefStretch(:)'; yDefStretch(:)']; % [ def1 x numOfRatios, def2 x numOfRatios, ...] x numOfSpeeds
            ratioStretched = ratioStretched(:)'; % repeat x numOfSpeedsByDefs
            speedStretched = [xSpeedStretched(:)'; ySpeedStretched(:)']; % [ speed1 x numOfRatiosByDefs, speed2 x numOfRatiosByDefs, ... ]
            driveParams = aol_drive_params(focalLength, optimalBaseFreq, xyDefStretched, ratioStretched, speedStretched, opWavelenVac);
        end
    end
    
end

