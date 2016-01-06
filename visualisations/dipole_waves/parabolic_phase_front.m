function create_phase_front()
tic
    function z = ConstantPhaseFunction(x)
        z = 0.2 * x.^2;
    end
    function dzdx = PhaseFunctionDerivative(x)
        dzdy = 0.4 * x;
    end

x0Test = 1;

kMod = 40;
x = -5:0.02:5;
z = -1:0.04:10;
[xMesh,zMesh] = meshgrid(x,z);
amp = Amplitude();
clf;
subplot(1,2,1);
MakeSubplot(abs(amp));
subplot(1,2,2);
MakeSubplot(real(amp));

    function MakeSubplot(ampFun)
        plot(-10:0.01:10, ConstantPhaseFunction(-10:0.01:10),'k','LineWidth', 2);
        hold on;
        surface(xMesh,zMesh,ampFun,'LineStyle', 'none');       
        plot([x0Test x0Test],[-10 10]);
        plot([-10 10], [25/x0Test, -25/x0Test]+(ConstantPhaseFunction(x0Test) + 5/2));
        plot([-10 10], [-4*x0Test, 4*x0Test]+(ConstantPhaseFunction(x0Test) - 2/5*x0Test*x0Test));
        theta = 0:0.1:2*pi;
        plot(0+2.5*cos(theta),2.5+2.5*sin(theta));
        axis equal;
        grid on;
        grid minor;
        axis([-5,5,-1,10]);
        hold off;
    end

    function amp = Amplitude()
        ampArray = integral(@Integrand, -1.4, 1.4, 'ArrayValued', true);
        %ampArray = Integrand(x0Test);
        amp = reshape(ampArray, size(xMesh));
    end

    function val = Integrand(x0)
        z0 = ConstantPhaseFunction(x0);
        theta = angle(1 + 1i * 2/5 * x0);
        k = [kMod .* -sin(theta); kMod .* cos(theta)];
        phase = k(1,:)*(xMesh(:)-x0) + k(2,:)*(zMesh(:)-z0);
        val = exp(1i*phase);
    end
toc
end

