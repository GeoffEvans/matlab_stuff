function [ f ] = GifOfPolarisation( amplitudeVector, titlePlot )
% Animates the given amplitude vector
% for left-hand polarised light, use amplitudeVector [1, exp(1i * pi/2)]
% for right-hand polarised light, use amplitudeVector [1, exp(1i * pi/2)]

T = 0:pi/20:2*pi;
Z = 0:pi/20:5*pi;
eFieldX = amplitudeVector(1) * exp(1i*Z);
eFieldY = amplitudeVector(2) * exp(1i*Z);

%Configure the graph
eFieldP = plot3(Z, eFieldX, eFieldY);
set(eFieldP,'YDataSource','eFieldX');
set(eFieldP,'ZDataSource','eFieldY');

xlabel('z position')
ylabel('x field')
zlabel('y field')
title(titlePlot ,'Color',[.6 0 0])
hold on
grid on
axis auto

f(1) = getframe(gcf);
[im,map] = rgb2ind(f(1).cdata, 256, 'nodither');

%Run over time
for j = 2:length(T)
    eFieldX = amplitudeVector(1) * exp(1i*(Z-T(j)));
    eFieldY = amplitudeVector(2) * exp(1i*(Z-T(j)));
    refreshdata(eFieldP,'caller');
    drawnow;
    f(j) = getframe(gcf);
    im(:,:,1,j) = rgb2ind(f(j).cdata, map, 'nodither');
end

hold off;
imwrite(im,map, strcat(titlePlot, '.gif'), 'DelayTime', 0, 'LoopCount', inf)

end

