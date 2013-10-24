function [ ] = PlotRays( rays )
    
startX = rays.starts(1,:);
startY = rays.starts(2,:);
endX = rays.stops(1,:);
endY = rays.stops(2,:);

line([startX; endX], [startY; endY], 'color', 'red');

end

