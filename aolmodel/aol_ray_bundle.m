classdef aol_ray_bundle < handle
    
    properties 
        numOfTimes
        numOfPositions
        numOfDrives
        numOfPerturbations
        numOfRays
        numOfRaysPerPerturbation
        
        t
        eff
        k
        xyz

        drives
        perturbsTheta
        perturbsPhi
    end
    
    methods
        function obj = aol_ray_bundle(microSecs,xyInputMm,drives,thetaPhiPerturbs)
            if nargin == 0
                microSecs = 0;
                xyInputMm = [0;0];
            end
            numOfAods = size(thetaPhiPerturbs{1},2);

            % plan of attack: put the perturbations in the third dimenson for easy mapping over rotations in the model.

            numOfPlanesZ = numOfAods*2+4; % input plane, 2 planes per AOD, focusPredicted, focusModel, exit plane
            obj.numOfTimes = length(microSecs);
            obj.numOfPositions = size(xyInputMm,2);
            obj.numOfPerturbations = size(thetaPhiPerturbs{1},1);
            obj.numOfDrives = length(drives);
            obj.numOfRays = obj.numOfTimes*obj.numOfPositions*obj.numOfDrives*obj.numOfPerturbations;
            obj.numOfRaysPerPerturbation = obj.numOfTimes*obj.numOfPositions*obj.numOfDrives;            
            
            repmatArray = [1,obj.numOfRaysPerPerturbation,obj.numOfPerturbations];
            obj.t = repmat(microSecs * 1e-6,repmatArray);
            obj.eff = zeros(numOfAods,obj.numOfRaysPerPerturbation,obj.numOfPerturbations);
            obj.k = repmat([0;0;1]*2*pi/aod3d.opWavelenVac,repmatArray); % input laser is orthogonal to AOD centre line
            obj.xyz = cell(numOfPlanesZ,1);
            xyzInitial = [xyInputMm*1e-3;zeros(1,size(xyInputMm,2))];
            obj.xyz{1} = obj.StretchPositionArray(xyzInitial);
            
            obj.drives = StretchDrives(drives);
            obj.perturbsTheta = thetaPhiPerturbs{1};
            obj.perturbsPhi = thetaPhiPerturbs{2};
        end
        
        function SetXyzNthAodFront(obj, n, xyzIn)
            obj.xyz{2*n} = xyzIn;
        end
        function SetXyzNthAodBack(obj, n, xyzIn)
            obj.xyz{2*n + 1} = xyzIn;
        end
        function xyzOut = GetXyzNthAodFront(obj, n)
            xyzOut = obj.xyz{2*n};
        end
        function xyzOut = GetXyzNthAodBack(obj, n)
            xyzOut = obj.xyz{2*n + 1};
        end
        function xyzOut = GetXyzLeavingAol(obj)
            xyzOut = obj.xyz{end - 3};
        end
    end
    methods (Access = private)
        function [ reshaped ] = StretchPositionArray(obj,array)
            reshaped = repmat(array,obj.numOfTimes,obj.numOfDrives);
            reshaped = reshape(reshaped,3,obj.numOfRaysPerPerturbation);
            reshaped = repmat(reshaped,[1,1,obj.numOfPerturbations]);
        end
        function [ reshaped ] = StretchDrives(obj,array)
            reshaped = repmat(array,obj.numOfTimes*obj.numOfPosititions);
            reshaped = reshape(reshaped,1,obj.numOfRaysPerPerturbation);
            reshaped = repmat(reshaped,[1,1,obj.numOfPerturbations]);
        end
    end
end

