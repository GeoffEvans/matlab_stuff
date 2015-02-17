classdef aol_drive_freqs < handle
    
    properties
        chirp
        baseFreq
    end
    
    methods
        function obj = aol_drive_freqs(baseFreq,chirp)
            % input size is numOfAods,numOfDrives
            if nargin == 0
                obj.chirp = 0;
                obj.baseFreq = 0;
            else
                numOfDrives = size(chirp,2);
                obj(numOfDrives) = aol_drive_freqs; % Preallocate object array
                for i = 1:numOfDrives
                    obj(i).chirp = chirp(:,i); % want columns for cat comma sep lists 
                    obj(i).baseFreq = baseFreq(:,i);
                end
            end
        end
    end
end

