function gaussian_waves()
% works in distance units of mm
% uses g defined by exp(ik/2 * x'gx)
% for simple astigmatisms, diag(g) is the inverse of q

minArea = 1e3;

wavelength = 800e-6;
vx = 0;
vy = 0;
l0 = 100;
l1 = 40;
l2 = 40;
l3 = 40;
l4 = 700;

f1 = 1/( (vx + 1)/(l1 + l2 + 2*l3 + 2*l4 + l1*vx + l2*vx) );
f2 = 1/( (vy + 1)/(l2 + l3 + 2*l4 + l2*vy + l3*vy) );
f3 = 1/( (1-vx)/(2*(l3 + l4)) );
f4 = 1/( (1-vy)/(2*l4) );

[zPoints, lenZPoints, g] = SetSizeOfSystem();

opticalElements = {struct('position',1)}; % input plane

AppendToOpticalElements( struct('position',find(zPoints == l0),'transform',@(gIn) CylindricalLen(gIn,l1+l2+l3+l4,pi)) );
AppendToOpticalElements( struct('position',find(zPoints == l0+l1),'transform',@(gIn) CylindricalLen(gIn,l2+l3+l4,pi/2)) );

%AppendToOpticalElements( struct('position',find(zPoints == l0),'transform',@(gIn) SphericalLens(gIn,l1+l2+l3+l4)) );

% AppendToOpticalElements( struct('position',find(zPoints == l0),'transform',@(gIn) CylindricalLen(gIn,f1,pi)) );
% AppendToOpticalElements( struct('position',find(zPoints == l0+l1),'transform',@(gIn) CylindricalLen(gIn,f2,pi/2)) );
% AppendToOpticalElements( struct('position',find(zPoints == l0+l1+l2),'transform',@(gIn) CylindricalLen(gIn,f3,0)) );
% AppendToOpticalElements( struct('position',find(zPoints == l0+l1+l2+l3),'transform',@(gIn) CylindricalLen(gIn,f4,-pi/2)) );
AppendToOpticalElements( struct('position',find(zPoints == l0+l1+l2+l3+l4),'transform',@(x) x) );
AppendToOpticalElements( struct('position',lenZPoints,'transform',@(x) x) );

close all;
hold on;

initialWidthX = 20;
initialWidthY = 20;
g{1} = Initial(initialWidthX,initialWidthY,inf,inf); % Input beam

for m = 1:length(opticalElements)-1
    prevElement = opticalElements{m};
    nextElement = opticalElements{m+1};
    CalculateGsBetweenElements(prevElement, nextElement)
    g{nextElement.position} = nextElement.transform(g{nextElement.position});
    
    fill3([1 1 -1 -1]*initialWidthX,[1 -1 -1 1]*initialWidthY,repmat(zPoints(prevElement.position),1,4),repmat(prevElement.position,1,4))
    alpha(0.2);
end
    
PlotGaussian();
hold off;
display(minArea);

    function AppendToOpticalElements(element)
        opticalElements = [opticalElements {element}];
    end

    function [zPoints, lenZPoints, g] = SetSizeOfSystem()
        zPoints = 0:10:1000; 
        lenZPoints = length(zPoints);
        g = cell(1,lenZPoints);
    end

    function CalculateGsBetweenElements(prevElement, nextElement)
        for n = prevElement.position+1:nextElement.position;
            g{n} = FreeSpace(g{prevElement.position},zPoints(n)-zPoints(prevElement.position));
        end
    end

    function PlotGaussian()
        thetaFromX = (-1:0.01:1)*pi;
        lenTheta = length(thetaFromX);
        waistAtZ = zeros(lenZPoints,lenTheta); % each row describes the waist at given z

        % find waist ellipse by ellipse down z
        for n2 = 1:lenZPoints
            G = g{n2};
            waistAtZ(n2,:) = FindWaistAtZ(G,thetaFromX);
        end 

        [thetaMat, zPointsMat] = meshgrid(thetaFromX, zPoints);
        s = surf(waistAtZ.*cos(thetaMat),waistAtZ.*sin(thetaMat),zPointsMat);
        set(s,'LineStyle','none');
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis([[-1,1]*initialWidthX,[-1,1]*initialWidthY,0,max(zPoints)])
        axis square
    end

    function gOut = Initial(waistX,waistY,rocX,rocY) 
        qModSqrX = 1/( (1/rocX).^2 + (wavelength/pi/waistX.^2).^2 );
        qX = (1/rocX - 1i * wavelength/pi/waistX.^2) * qModSqrX;
        qModSqrY = 1/( (1/rocY).^2 + (wavelength/pi/waistY.^2).^2 );
        qY = (1/rocY - 1i * wavelength/pi/waistY.^2) * qModSqrY;
        q = [qX,0;0,qY];
        gOut = inv(q);
    end

    function gOut = FreeSpace(gIn,opticalPathLength)
        gOut = inv( inv(gIn) + opticalPathLength*eye(2) );
    end

    function gOut = CylindricalLen(gIn,focalLength,angleFromX)
        c = cos(angleFromX);
        s = sin(angleFromX);
        f = -[c*c,c*s;c*s,s*s]/focalLength;
        gOut = gIn + f;
    end

    function gOut = SphericalLens(gIn,focalLength)
        gOut = CylindricalLen(gIn,focalLength,0);
        gOut = CylindricalLen(gOut,focalLength,pi/2);
    end
    
    function [rMajor, rMinor, angleMajor] = Waist(g)
        waistMatrix = imag(g);
        [vectorsInCols,valuesMat] = eigs(waistMatrix);
        [vectorsInCols,values] = sortem(vectorsInCols,valuesMat);
        angleMajor = angle(vectorsInCols(1,2) + 1i * vectorsInCols(2,2));
        rMajor = sqrt(wavelength/pi/values(2));
        rMinor = sqrt(wavelength/pi/values(1));
        area = rMajor*rMinor*pi;
        if minArea > area
            minArea = area;
        end
    end

    function waist = FindWaistAtZ(g,theta)
        [rMajor, rMinor, angleMajor] = Waist(g);
        waist = (rMajor*rMinor + theta*0)./sqrt( (rMinor*cos(theta-angleMajor)).^2 + (rMajor*sin(theta-angleMajor)).^2 );
    end
    
	function phaseFront = PhaseFront(q)
        % do something clever to plot phase fronts inside waist
    end
    
end