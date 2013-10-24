function [ ] = PlotInterfaces( interfaces )

global maxY minY stepY;

Y = minY:stepY:maxY;
for j = 1:length(interfaces) 
    i = interfaces(j);
    plot(i.shape(Y) + i.position, Y, 'color', 'black');
end
global opticalLength;
axis([0 opticalLength minY maxY])


end

