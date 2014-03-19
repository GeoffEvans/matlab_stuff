classdef aol_drive < handle
    
    properties
        chirp
        baseFreq
    end
    
    methods
        function obj = aol_drive(baseFreq,chirp)
            if nargin == 0
                obj.chirp = 0;
                obj.baseFreq = 0;
            else
                numOfDrives = size(chirp,1);
                obj(numOfDrives) = aol_drive; % Preallocate object array
                for i = 1:numOfDrives
                    obj(i).chirp = chirp(i,:);
                    obj(i).baseFreq = baseFreq(i,:);
                end
            end
        end
    end
end

