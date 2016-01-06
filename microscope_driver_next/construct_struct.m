function s = construct_struct(varargin)
    s = struct();
    for n = 1:nargin
        s.(inputname(n)) = varargin{n};
    end
end
