function plot_acoustic_indicatrix()
% For various wavevector directions, l, find the phase velocity eigenvalues
% and displacement velocity eigenvectors.

% Orientations of input plane waves
thetaRange = 0:pi/90:pi;
phiRange = 0:pi/90:2*pi;
[thetaMesh, phiMesh] = meshgrid(thetaRange, phiRange);
[ unorderedVels, unorderedDispVecs ] = teo2.find_v_ac(thetaMesh, phiMesh);

[ velocities, displacementVectors ] =...
    cellfun(@SortModesByVelocity, unorderedDispVecs, unorderedVels, 'UniformOutput', false);

n1Mesh = ConvertVelocityToIndex(1);
n2Mesh = ConvertVelocityToIndex(2);
n3Mesh = ConvertVelocityToIndex(3);

[X1,Y1,Z1] = ConvertPolarToCartesian(n1Mesh);
[X2,Y2,Z2] = ConvertPolarToCartesian(n2Mesh);
[X3,Y3,Z3] = ConvertPolarToCartesian(n3Mesh);
[c1,c2,c3] = AssignColoursToHowTransverse();

PlotIndicatrixSurfaces();

    function mesh = ConvertVelocityToIndex(i) 
        mesh = cellfun(@(v) 1./v(i), velocities);
    end

    function [descendingArrayOfVels, sortedVectors] = SortModesByVelocity(dispVectors, vels) 
        [descendingArrayOfVels, indicesOfSortedArray] = sort(vels,'descend');
        sortedVectors = dispVectors(:,indicesOfSortedArray); % arrange the columns in this order
    end

    function [c1,c2,c3] = AssignColoursToHowTransverse()
        HowTransverse1 = @(t,p,v) HowTransverse(t,p,v,1);
        HowTransverse2 = @(t,p,v) HowTransverse(t,p,v,2);
        HowTransverse3 = @(t,p,v) HowTransverse(t,p,v,3);
        c1 = cellfun(HowTransverse1, num2cell(thetaMesh), num2cell(phiMesh), displacementVectors);
        c2 = cellfun(HowTransverse2, num2cell(thetaMesh), num2cell(phiMesh), displacementVectors);
        c3 = cellfun(HowTransverse3, num2cell(thetaMesh), num2cell(phiMesh), displacementVectors);
        
        function f = HowTransverse(theta,phi,threeDispVecs,i) 
            unitNormal = [sin(theta).*cos(phi) sin(theta).*sin(phi) cos(theta)];
            dispVec = threeDispVecs(:,i);
            f = 1 - dot(unitNormal, dispVec).^2;
        end
    end

    function PlotIndicatrixSurfaces()
        s1 = surf(X1,Y1,Z1,c1,'EdgeColor','none');
        hold on;
        s2 = surf(X2,Y2,Z2,c2,'EdgeColor','none');
        s3 = surf(X3,Y3,Z3,c3,'EdgeColor','none');
        % set(gcf, 'renderer', 'zbuffer'); % Needed for colorbar to work on ucbtgje
        % % but breaks alpha.
        grid on;
        axis square;
        xlabel('x');
        ylabel('y');
        zlabel('z optic axis');
        alpha(s1, 0.2);
        alpha(s2, 0.2);
        alpha(s3, 0.2);
        alphamap('rampdown');
        camlight right;
        lighting phong;
    end

    function [X,Y,Z] = ConvertPolarToCartesian(radiusMesh)
        X = radiusMesh .* sin(thetaMesh).* cos(phiMesh);
        Y = radiusMesh .* sin(thetaMesh).* sin(phiMesh);
        Z = radiusMesh .* cos(thetaMesh);
    end

end