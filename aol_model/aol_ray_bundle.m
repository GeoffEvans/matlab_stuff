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
        
        xyFocusModel
        zFocusModel
        zFocusPredicted
        
        aodCentres
        drives
        perturbations
    end
    
    methods
        function obj = aol_ray_bundle(microSecs,xyInputMm,drives,aodCentres,zFocusPredicted,aolPerturbs)
            if nargin == 0
                microSecs = 0;
                xyInputMm = [0;0];
            end
            obj.numOfAods = aolPerturbs.numOfAods;
            
            % plan of attack: put the perturbations in the third dimension of xyz,k,eff for easy mapping over rotations in the model.
            
            numOfPlanesZ = obj.numOfAods*2+4; % input plane, 2 planes per AOD, focusPredicted, focusModel, exit plane
            obj.numOfTimes = length(microSecs);
            obj.numOfPositions = size(xyInputMm,2);
            obj.numOfPerturbations = aolPerturbs.numOfPerturbations;
            obj.numOfDrives = length(drives);
            obj.numOfRays = obj.numOfTimes*obj.numOfPositions*obj.numOfDrives*obj.numOfPerturbations;
            obj.numOfRaysPerPerturbation = obj.numOfTimes*obj.numOfPositions*obj.numOfDrives;
                      
            obj.eff = zeros(obj.numOfAods,obj.numOfRaysPerPerturbation,obj.numOfPerturbations);
            
            repmatArrayK = [1,obj.numOfRaysPerPerturbation,obj.numOfPerturbations];
            obj.k = repmat([0;0;1]*2*pi/aod3d.opWavelenVac,repmatArrayK); % input laser is orthogonal to AOD centre line
            
            repmatArrayT = [1,obj.numOfPositions*obj.numOfDrives,obj.numOfPerturbations];
            obj.t = repmat(microSecs * 1e-6,repmatArrayT);
            
            obj.xyz = cell(numOfPlanesZ,1);
            xyzInitial = [xyInputMm*1e-3;zeros(1,size(xyInputMm,2))];
            obj.xyz{1} = obj.StretchPositionArray(xyzInitial);
            
            obj.zFocusPredicted = zFocusPredicted;
            obj.drives = obj.StretchDrives(drives);
            obj.aodCentres = cellfun(@obj.StretchDrives, aodCentres, 'UniformOutput', 0); % one centre for each drive (pairDefRat and xyDef)
            obj.perturbations = aolPerturbs;
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
            drivesLocal = obj.drives;
            bf = [drivesLocal.baseFreq];
            bf = reshape(bf,[obj.numOfAods,obj.numOfRaysPerPerturbation,obj.numOfPerturbations]);
            bf = bf(nthAod,:,:);
        end
        
        function c = ChirpForNthAod(obj, nthAod)
            drivesLocal = obj.drives;
            c = [drivesLocal.chirp];
            c = reshape(c,[obj.numOfAods,obj.numOfRaysPerPerturbation,obj.numOfPerturbations]);
            c = c(nthAod,:,:);
        end
        
        function vectorsOut = ApplyPerturbationMatricesToVectors(obj, MapPerturbationToMatrix, vectorsIn, nthAod)
            vectorsOut = zeros(size(vectorsIn));
            [thetaPerturbs,phiPerturbs] = obj.perturbations.GetPerturbationsForAod(nthAod);
            for m = 1:obj.numOfPerturbations
                vectorsOut(:,:,m) = MapPerturbationToMatrix(nthAod,thetaPerturbs(m),phiPerturbs(m)) * vectorsIn(:,:,m);
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

