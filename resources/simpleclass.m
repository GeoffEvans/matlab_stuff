classdef simpleclass
    
    properties
        propA 
        propB
    end
    
    methods
%         function obj = simpleclass()
%             obj.propA = 1;
%         end
        function obj = setb(obj,v)
            obj.propB = v;
        end
    end
end

