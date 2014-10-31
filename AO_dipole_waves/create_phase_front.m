function create_phase_front()
tic
    a0 = 0;
    a1 = 0;
    a2 = 0.2;
    a3 = 0.0;
    function z = ConstantPhaseFunction(x)
        z = a3*x.^3 + a2*x.^2 + a1*x + a0;
    end
    function dzdx = PhaseFunctionDerivative(x)
        dzdx = 0*a3*x.^2 + 2*a2*x + 0*a1;
    end

x0test = 1;

kMod = 40;
x = -5:0.05:5;
z = -1:0.05:10;
[xMesh,zMesh] = meshgrid(x,z);
frontRange = -1.4:0.01:1.4;
amp = Amplitude();
clf;
subplot(1,2,1);
MakeSubplot(abs(amp));
subplot(1,2,2);
MakeSubplot(real(amp));

    function MakeSubplot(ampFun)
        plot(frontRange, ConstantPhaseFunction(frontRange),'k','LineWidth', 2);
        hold on;
        surface(xMesh,zMesh,ampFun,'LineStyle', 'none');       
        PlotGuideLines()
        axis equal;
        grid on;
        grid minor;
        axis([min(x),max(x),min(z),max(z)]);
        hold off;
        
        function PlotGuideLines()
            z0test = ConstantPhaseFunction(x0test);
            deriv = PhaseFunctionDerivative(x0test);
            plot([-10 10], [(x0test + 10)./deriv, (x0test - 10)./deriv] + z0test);
            plot([-10 10], [(-10 - x0test).*deriv, (10 - x0test).*deriv] + z0test);
            if a1 == 0 && a0 == 0 && a3 == 0
                theta = 0:0.1:2*pi;
                plot(0+0.5/a2*cos(theta),0.5/a2*(1+sin(theta)));
            end
        end
    end

    function amp = Amplitude()
        ampArray = integral(@PlaneWaveAlongFront, min(frontRange), max(frontRange), 'ArrayValued', true);
        amp = reshape(ampArray, size(xMesh));
    end

    function val = PlaneWaveAlongFront(x0)
        z0 = ConstantPhaseFunction(x0);
        theta = angle(1 + 1i * PhaseFunctionDerivative(x0));
        k = [kMod .* -sin(theta); kMod .* cos(theta)];
        phase = k(1,:)*(xMesh(:)-x0) + k(2,:)*(zMesh(:)-z0);
        val = exp(1i*phase);
    end
toc
end

