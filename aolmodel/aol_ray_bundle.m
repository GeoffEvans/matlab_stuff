classdef aol_ray_bundle < handle
% Principally a class to manage information on each ray ->
% time,position,efficiency,momentum,centreOfAod,driveFreqs,perturbationOfAod

    properties
        numOfAods
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
        
        zFocusModel
        zFocusPredicted
        
        aodCentres
        drives
        
        perturbsTheta
        perturbsPhi
    end
    
    methods
        function obj = aol_ray_bundle(microSecs,xyInputMm,drives,aodCentres,zFocusPredicted,thetaPhiPerturbs)
            if nargin == 0
                microSecs = 0;
                xyInputMm = [0;0];
            end
            obj.numOfAods = size(thetaPhiPerturbs{1},2);
            
            % plan of attack: put the perturbations in the third dimenson for easy mapping over rotations in the model.
            
            numOfPlanesZ = obj.numOfAods*2+4; % input plane, 2 planes per AOD, focusPredicted, focusModel, exit plane
            obj.numOfTimes = length(microSecs);
            obj.numOfPositions = size(xyInputMm,2);
            obj.numOfPerturbations = size(thetaPhiPerturbs{1},1);
            obj.numOfDrives = length(drives);
            obj.numOfRays = obj.numOfTimes*obj.numOfPositions*obj.numOfDrives*obj.numOfPerturbations;
            obj.numOfRaysPerPerturbation = obj.numOfTimes*obj.numOfPositions*obj.numOfDrives;
                      
            obj.eff = zeros(obj.numOfAods,obj.numOfRaysPerPerturbation,obj.numOfPerturbations);
            repmatArray = [1,obj.numOfRaysPerPerturbation,obj.numOfPerturbations];
            obj.k = repmat([0;0;1]*2*pi/aod3d.opWavelenVac,repmatArray); % input laser is orthogonal to AOD centre line
            repmatArray = [1,obj.numOfPositions*obj.numOfDrives,obj.numOfPerturbations];
            obj.t = repmat(microSecs * 1e-6,repmatArray);
            
            obj.xyz = cell(numOfPlanesZ,1);
            xyzInitial = [xyInputMm*1e-3;zeros(1,size(xyInputMm,2))];
            obj.xyz{1} = obj.StretchPositionArray(xyzInitial);
            
            obj.zFocusPredicted = zFocusPredicted;
            obj.aodCentres = aodCentres;
            obj.drives = obj.StretchDrives(drives);
            obj.perturbsTheta = thetaPhiPerturbs{1};
            obj.perturbsPhi = thetaPhiPerturbs{2};
        end
        
        function SetXyzNthAodFront(obj, n, xyzIn)
            obj.xyz{2*n} = xyzIn;
        end
        function SetXyzNthAodBack(obj, n, displacementInCrystal)
            obj.xyz{2*n + 1} = obj.xyz{2*n} + displacementInCrystal;
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
        
        function bf = BaseFreqForNthAod(obj, nthAod)
            % TODO
        end
        
        function c = ChirpForNthAod(obj, nthAod)
            % TODO 
        end
        
        function vectorsOut = ApplyPerturbationMatricesToVectors(obj, MapPerturbationToMatrix, vectorsIn, nthAod)
            vectorsOut = zeros(size(vectorsIn));
            for m = 1:obj.numOfPerturbations
                phiAod = obj.perturbsPhi(m,nthAod);
                thetaAod = obj.perturbsPhi(m,nthAod);
                vectorsOut(:,:,m) = MapPerturbationToMatrix(nthAod,thetaAod,phiAod) * vectorsIn(:,:,m);
            end
        end
    end
    methods (Access = private)
        function [ reshaped ] = StretchPositionArray(obj,array)
            reshaped = stretch(array,obj.numOfTimes);
            reshaped = repmat(reshaped,[1,obj.numOfDrives,obj.numOfPerturbations]);
        end
        function [ reshaped ] = StretchDrives(obj,array)
            reshaped = stretch(array,obj.numOfTimes*obj.numOfPositions);
            reshaped = repmat(reshaped,[1,1,obj.numOfPerturbations]);
        end
    end
end

