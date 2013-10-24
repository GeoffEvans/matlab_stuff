% This script animates the effect of a quarter wave plate at different
% angles to the major axis of the elliptical polarisation.

% Define ellipse params
xMax = 5;
yMax = 2;

S = 0:0.1:2*pi; % overall phase
X = real(xMax * exp(1i * S));
Y = real(yMax * exp(1i * (S + pi/2))); % so this is left-handed in theory

O = 0; % angleToMajorAxis
[X2,Y2,X3L,Y3L,X4L,Y4L,X3R,Y3R,X4R,Y4R,pL,pR] = CalculateCoords(xMax,yMax,O,S);
[majorAxis,~,~] = RotationCoords(O,X,Y);
fig = figure();
set(gcf, 'Position', [50 50 800 500]);

subplot(2,3,1); % View ellipse along principal axes
PlotEllipse(X,Y);
xlabel('E major');
ylabel('E minor');
line([-5 5],[0 0], 'color', 'black');
title('Handedness: red-right, blue-left')

subplot(2,3,4); % View ellipse aligned against fast-slow axes
p2 = PlotEllipse(X2,Y2);
xlabel('E fast');
ylabel('E slow');
line([-5 5],[0 0], 'color', 'black');
l = line(majorAxis(1,:), majorAxis(2,:), 'color', 'black');

subplot(2,3,5); % View NEW LHP ellipse aligned against fast-slow axes
p3 = PlotEllipse(X3L,Y3L);
xlabel('E fast')
ylabel('E slow')
subplot(2,3,2); % View NEW LHP ellipse along principal axes
p4 = PlotEllipse(X4L,Y4L);
xlabel('E major');
ylabel('E minor');
title('Initially left-hand');
subplot(2,3,6); % View NEW RHP ellipse aligned against fast-slow axes
p5 = PlotEllipse(X3R,Y3R);
xlabel('E fast')
ylabel('E slow')
subplot(2,3,3); % View NEW RHP ellipse along principal axes
p6 = PlotEllipse(X4R,Y4R);
xlabel('E major');
ylabel('E minor');
title('Initially right-handed');

j = 1;
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata, 256, 'nodither');

for O = 0:pi/90:pi;
    [X2,Y2,X3L,Y3L,X4L,Y4L,X3R,Y3R,X4R,Y4R,lPrh,rPrh] = CalculateCoords(xMax,yMax,O,S);
    [majorAxis, Xe, Ye] = RotationCoords(O,X,Y);
    refreshdata(p2,'caller');
    refreshdata(p3,'caller');
    refreshdata(p4,'caller');    
    refreshdata(p5,'caller');
    refreshdata(p6,'caller');
    if (lPrh) % Left-handed ellipse is now polarised right-handed (prh)
        set(p3,'color','red');
        set(p4,'color','red');
    else
        set(p3,'color','blue');
        set(p4,'color','blue');
    end
    if (rPrh)
        set(p5,'color','red');
        set(p6,'color','red');
    else
        set(p5,'color','blue');
        set(p6,'color','blue');
    end
    set(l,'XData',majorAxis(1,:));
    set(l,'YData',majorAxis(2,:));
    drawnow;
    j = j+1;
    f = getframe(gcf);
    im(:,:,1,j) = rgb2ind(f.cdata, map, 'nodither');
end

imwrite(im,map, strcat('quarter-ellipse', '.gif'), 'DelayTime', 0, 'LoopCount', inf)