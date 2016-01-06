function [ m ] = mag( v )
% Maps magnitude over a row of vectors
m = dot(v,v);
m = sqrt(m);
end


