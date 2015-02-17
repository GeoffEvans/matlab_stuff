classdef Person
    %PERSON Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
    end
    
    methods
        function obj = Person(name)
            obj.name = obj.change_name(name);
        end
        function give_name(blah)
            display(blah.name)
        end
        function up = change_name(blah, name)
            up = upper(name);
        end
    end
    
end

