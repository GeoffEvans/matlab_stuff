function chirped()

w_op = 5;
k_op = 20;

V_ac = 1;
w_ac0 = 6;
w_ac1 = 0;
w_ac_period = 10;

figure();%figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w'); % set figure background to white
p = MakePlot();

AnimatePlots();

    function p = MakePlot()
        [xMesh,zMesh,imageMesh] = Amplitudes(0);
        p = MakePlot(imageMesh,'+1');
                
        function p = MakePlot(colourMesh,titleStr)
            p = pcolor(xMesh,zMesh,colourMesh);
            set(p,'EdgeColor', 'none');
            xlabel('[110]')
            ylabel('[001]')
            title(titleStr)
        end
    end

    function AnimatePlots()
        slides = 300;
        %f(1) = getframe(gcf);
        %[~,map] = rgb2ind(f(1).cdata, 256, 'nodither');
        for k = 1:slides
            UpdatePlots(k);
            %f(k) = getframe(gcf);
            %im(:,:,1,k) = rgb2ind(f(k).cdata, map, 'nodither');
        end
        %imwrite(im,map, strcat('AOD scattering', '.gif'), 'DelayTime', 0, 'LoopCount', inf)

        function UpdatePlots(n)
            t = 2 * pi / w_op * n / slides * 20; % last number gives number of optical periods
            [~,~,imageMesh] = Amplitudes(t);
            set(p,'CData', imageMesh);
            drawnow;
        end
    end

    function [xMesh,zMesh,imageMesh] = Amplitudes(t)
        % This would be much nicer done with a piecewise function...
        X = -10:0.01:5;
        Z = -1:0.01:0.5;
        
        [xMesh,zMesh] = meshgrid(X,Z);
        imageMesh = ImageAmp(xMesh,zMesh,t);
    end

    function amp = ImageAmp(x,z,t)
        amp = PlusOneAmp(x,z,t); %initialise
        opticFilter = z < -0.5;
        acousticFilter = x < 5;
        amp(opticFilter) = OpticAmp(x(opticFilter),z(opticFilter),t);
        amp(acousticFilter) = AcousticAmp(x(acousticFilter),z(acousticFilter),t);
        amp(acousticFilter & opticFilter) = 0;
    end
    function amp = PlusOneAmp(x,z,t)
        amp = cos( (w_op + AcousticFreq(t,x)).*t - AcousticFreq(t,x)/V_ac.*x - k_op.*z );
    end
    function amp = OpticAmp(x,z,t)
        amp = cos( w_op.*t - k_op.*z );
    end
    function w = AcousticFreq(t,x)
        w = w_ac0 + w_ac1*sawtooth((t - x/V_ac)/w_ac_period);
    end
    function amp = AcousticPlaneWave(kz,x,z,t)
        k = AcousticFreq(t,x)/V_ac;
        kx = sqrt(k.^2 - kz.^2);
        amp = cos( AcousticFreq(t,x).*t - kx.*x - kz.*z );
    end
    function amp = AcousticAmp(x,z,t)
        k = AcousticFreq(t,x)/V_ac;
        waves = 0*x;
        for m = 0:10
            kz = m*k/11 - k/2;
            waves = waves + 0.1*AcousticPlaneWave(kz,x,z,t).*sinc(kz*3/pi/2);
        end
        amp = waves;
    end
end

