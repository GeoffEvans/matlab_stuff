function plot_optical_indicatrix()
% For various wavevector directions, find the refractive indices and their respective polarisations

% Orientations of input plane waves
thetaRadians = -pi/2:pi/360:pi/2;
thetaDegrees = thetaRadians*180/pi;
[ n_ord, n_ext, p_ord, p_ext ] = teo2.find_n_op(thetaRadians);

PlotIndexAgainstTheta();
PlotSurfaceEllipsoids();

    function PlotIndexAgainstTheta()
        subplot(1, 2, 1);
        hold on;
        grid on;
        grid minor;
        xlabel('theta');
        ylabel('n');
        PlotApproximateEllipsoids();
        plot(thetaDegrees,real(n_ord), 'color', 'red');
        plot(thetaDegrees,real(n_ext));
        hold off;
    end

    function PlotApproximateEllipsoids()
        delta = n_ord(1).^2 * teo2.activityVector / 2; %(1.83 Xu & St)
        ordNumerator = ( cos(thetaRadians)/((1 - delta)*n_ord(1)) ).^2 + ( sin(thetaRadians)/n_ord(1) ).^2;
        extNumerator = ( cos(thetaRadians)/((1 + delta)*n_ord(1)) ).^2 + ( sin(thetaRadians)/n_ext(1) ).^2;
        nOrdApprox = ordNumerator .^ -0.5;
        nExtApprox = extNumerator .^ -0.5;
        plot(thetaDegrees,nOrdApprox, 'color', 'green');
        plot(thetaDegrees,nExtApprox, 'color', 'magenta');
    end

    function f = AssignColourToPolarisation(pol)
        % Calculate the handedness of the polarisations: negative phase implied
        % left-handed, very big or small ratio implies linear
        RightHanded = @(pol) angle(pol) > 0;
        f = (RightHanded(pol)-0.5) .* log(4 + log(4 + abs(pol)));
    end


    function PlotSurfaceEllipsoids()
        phiRange = 0:pi/360:2*pi;
        phiLength = length(phiRange);
        [thetaMesh, phiMesh] = meshgrid(thetaRadians, phiRange);
        n1Mesh = ones(phiLength,1) * n_ord;
        n2Mesh = ones(phiLength,1) * n_ext;
        [X1,Y1,Z1] = ConvertPolarToCartesian(n1Mesh);
        [X2,Y2,Z2] = ConvertPolarToCartesian(n2Mesh);
        c1Mesh = ones(phiLength,1) * AssignColourToPolarisation(p_ord);
        c2Mesh = ones(phiLength,1) * AssignColourToPolarisation(p_ext);

        subplot(1, 2, 2);
        s = surf(X1,Y1,Z1,c1Mesh,'EdgeColor','none');
        alpha(s, 0.3);
        hold on;
        t = surf(X2,Y2,Z2,c2Mesh,'EdgeColor','none');
        alpha(t, 0.2);
        grid on;
        axis square;
        xlabel('x');
        ylabel('y');
        zlabel('z optic axis');
                
        function [X,Y,Z] = ConvertPolarToCartesian(radiusMesh)
            X = radiusMesh .* sin(thetaMesh).* cos(phiMesh);
            Y = radiusMesh .* sin(thetaMesh).* sin(phiMesh);
            Z = radiusMesh .* cos(thetaMesh);
        end
    end
end

