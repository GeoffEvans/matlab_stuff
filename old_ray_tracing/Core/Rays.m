classdef Rays
    
    properties
        starts
        gradients
        stops
    end
    
    methods
        function r = Rays(starts, gradients, stops)
            r.starts = starts;
            r.gradients = gradients;
            r.stops = stops;
        end
    end
    
end

