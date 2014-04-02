classdef aol_perturbations < handle
    
    properties
        theta
        phi
        numOfPerturbations = 1;
        numOfAods
    end
    
    methods
        function obj = aol_perturbations(theta,phi)
            if length(theta) ~= length(phi)
                error('input sizes must be equal');
            end
            obj.numOfAods = length(theta);
            obj.theta = theta(:);
            obj.phi = phi(:);
        end
        function AddPerturbation(obj, theta, phi)
            if length(theta) ~= length(phi)
                error('input sizes must be equal');
            end
            if obj.numOfAods ~= length(theta)
                error('number of AODs mismatch');
            end
            
            obj.theta = [obj.theta, theta(:)];
            obj.phi = [obj.phi, phi(:)];
            obj.numOfPerturbations = obj.numOfPerturbations + 1;
        end
        function [thetaN,phiN] = GetPerturbationsForAod(obj,n)
            thetaN = obj.theta(n,:);
            phiN = obj.phi(n,:);
        end
    end
end

