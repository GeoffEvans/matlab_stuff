classdef Interface
    
    properties
        position
        shape
        refractiveIndex
    end
    
    methods
      function i = Interface(position, shape, refractiveIndex)
         i.position = position;
         i.shape = shape;
         i.refractiveIndex = refractiveIndex;
      end
    end
    
end

