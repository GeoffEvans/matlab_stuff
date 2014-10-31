function gaussian_waves()
% works in distance units of mm
% uses g defined by exp(ik/2 * x'gx)
% for simple astigmatisms, diag(g) is the inverse of q

% Remember this function uses the principal that
% 1) The Gaussian phase term can be wrapped up as (x y)(A B;C D)(x;y) and so
% passage through a lens is trivial to add to it
% 2) The Gaussian phase term propagating through free space goes, for a simple astigmatic
% beam, as (A B;C D) -> (A B;C D) + zI. A general astigmatic beam can be written as a
% (complex) rotation, R,  of a simple astigmatic beam and since rotations are
% self-inverse R(A B;C D)R' -> R(A B;C D)R' + zI.

minArea = 1e3;
wavelength = 800e-6;
initialWidthX = 20;
initialWidthY = 10;

[zPoints, g] = SetSizeOfSystem();
opticalElements = SpecifyOpticalElements(zPoints);

close all;
hold on;
PropagateGaussian(opticalElements);
PropagateRays(opticalElements);
hold off;

    function PropagateGaussian(opticalElements)
        g{1} = Initial(initialWidthX,initialWidthY,inf,inf); % Input beam
        
        for m = 1:length(opticalElements)-1
            prevElement = opticalElements(m);
            nextElement = opticalElements(m+1);
            CalculateGsBetweenElements(prevElement, nextElement)
            g{nextElement.position} = nextElement.transform(g{nextElement.position});
            
            fill3([1 1 -1 -1]*initialWidthX,[1 -1 -1 1]*initialWidthY,repmat(zPoints(prevElement.position),1,4),repmat(prevElement.position,1,4))
            alpha(0.2);
        end
        
        PlotGaussian();
        display(minArea);
        
        function gOut = Initial(waistX,waistY,rocX,rocY)
            qModSqrX = 1/( (1/rocX).^2 + (wavelength/pi/waistX.^2).^2 );
            qX = (1/rocX - 1i * wavelength/pi/waistX.^2) * qModSqrX;
            qModSqrY = 1/( (1/rocY).^2 + (wavelength/pi/waistY.^2).^2 );
            qY = (1/rocY - 1i * wavelength/pi/waistY.^2) * qModSqrY;
            q = [qX,0;0,qY];
            gOut = inv(q);
        end
        
        function CalculateGsBetweenElements(prevElement, nextElement)
            for n = prevElement.position+1:nextElement.position;
                g{n} = FreeSpace(g{prevElement.position},zPoints(n)-zPoints(prevElement.position));
            end
            
            function gOut = FreeSpace(gIn,opticalPathLength)
                gOut = inv( inv(gIn) + opticalPathLength*eye(2) );
            end
        end
        
        function waist = FindWaistAtZ(g,theta)
            [rMajor, rMinor, angleMajor] = Waist(g);
            waist = (rMajor*rMinor + theta*0)./sqrt( (rMinor*cos(theta-angleMajor)).^2 + (rMajor*sin(theta-angleMajor)).^2 );
            
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
        end
        
        function PlotGaussian()
            thetaFromX = (-1:0.01:1)*pi;
            thetaLength = length(thetaFromX);
            waistAtZ = zeros(length(zPoints),thetaLength); % each row describes the waist at given z
            
            % find waist ellipse by ellipse down z
            for n2 = 1:length(zPoints)
                G = g{n2};
                waistAtZ(n2,:) = FindWaistAtZ(G,thetaFromX);
            end
            
            [thetaMat, zPointsMat] = meshgrid(thetaFromX, zPoints);
            s = surf(waistAtZ.*cos(thetaMat),waistAtZ.*sin(thetaMat),zPointsMat);
            set(s,'LineStyle','none');
            alpha(0.3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis([[-1,1]*initialWidthX,[-1,1]*initialWidthY,0,max(zPoints)])
            axis square
        end
    end
    function PropagateRays(opticalElements)
        numOpticElems = length(opticalElements);
        rayCount = 10;
        rayPositions = zeros(3,rayCount,numOpticElems); % each matrix corresponds to an optical element
        rayDirections = zeros(2,rayCount); % only x,y directions
        
        rayPositions(:,:,1) = DefineInitialRays();
        
        for n = 2:numOpticElems
            rayPositions(:,:,n) = PropagateRayFreeSpace(n,rayPositions,rayDirections);
            rayDirections = DeflectAtNthElement(n,rayPositions,rayDirections);          
        end
        
        function rays = DefineInitialRays()
            thetas = linspace(-pi,pi,rayCount);
            rays = [initialWidthX * cos(thetas); initialWidthY * sin(thetas); ones(1,rayCount)];
        end 
        
        function newRayPositions = PropagateRayFreeSpace(n,rayPositions,rayDirections)
            zDisplacement = LocationOfNthElement(n) - LocationOfNthElement(n-1);
            displacement = [rayDirections * zDisplacement; repmat(zDisplacement,1,rayCount)]; % note in air, rayDirections has mag 1
            newRayPositions = rayPositions(:,:,n-1) + displacement; 

            x = [newRayPositions(1,:) ; rayPositions(1,:,n-1)];
            y = [newRayPositions(2,:) ; rayPositions(2,:,n-1)];
            z = [newRayPositions(3,:) ; rayPositions(3,:,n-1)];
            line(x,y,z);
            
            function loc = LocationOfNthElement(n)
                loc = zPoints(opticalElements(n).position);
            end
        end
        
        function newRayDirections = DeflectAtNthElement(n,rayPositions,rayDirections)
            product = zeros(2,rayCount);
            for m = 1:rayCount
                product(:,m) = opticalElements(n).transform(zeros(2,2)) * rayPositions(1:2,m,n);
            end
            newRayDirections = rayDirections + product;
        end
    end
end

function [zPoints, g] = SetSizeOfSystem()

zPoints = 0:10:2500;
g = cell(1,length(zPoints));

end

function opticalElems = SpecifyOpticalElements(zPoints)
l0 = 100;
l1 = 40;
l2 = 40;
l3 = 40;
l4 = 700;

opticalElems = struct('position',1,'transform',@(x) x); % input plane

AppendToOpticalElements( struct('position',find(zPoints == l0),'transform',@(gIn) CylindricalLen(gIn,l1+l2+l3+2*l4,pi)) );
AppendToOpticalElements( struct('position',find(zPoints == l0+l4),'transform',@(gIn) CylindricalLen(gIn,l2+l3+2*l4,pi/2.3)) );

AppendToOpticalElements( struct('position',length(zPoints),'transform',@(x) x) ); % end plane

    function AppendToOpticalElements(element)
        opticalElems = [opticalElems element];
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
end