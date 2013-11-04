function monochromatic()

w_op = 5;
w_ac = 0.1;
k_op = 20;
k_ac = 1;

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w'); % set figure background to white
[p1,p2,p3] = MakePlots();

AnimatePlots();

    function [p1,p2,p3] = MakePlots()
        [xMesh,zMesh,totalProductMesh,totalPlusOneMesh,totalMinusOneMesh] = Amplitudes(0);
        p1 = MakeSubplot(1,totalProductMesh,'superposition');
        p2 = MakeSubplot(2,totalPlusOneMesh,'+1');
        p3 = MakeSubplot(3,totalMinusOneMesh,'-1');
        
        function p = MakeSubplot(i,colourMesh,titleStr)
            subplot(3,1,i);
            p = pcolor(xMesh,zMesh,colourMesh);
            set(p,'EdgeColor', 'none');
            xlabel('[110]')
            ylabel('[001]')
            title(titleStr)
        end
    end

    function AnimatePlots()
        
        slides = 800;
        %f(1) = getframe(gcf);
        %[~,map] = rgb2ind(f(1).cdata, 256, 'nodither');
        for k = 1:slides
            UpdatePlots(k);
            %f(k) = getframe(gcf);
            %im(:,:,1,k) = rgb2ind(f(k).cdata, map, 'nodither');
        end
        %imwrite(im,map, strcat('AOD scattering', '.gif'), 'DelayTime', 0, 'LoopCount', inf)

        function UpdatePlots(k)
            t = 2 * pi / w_op * k / slides * 20; % last number gives number of optical periods
            [~,~,totalProductMesh,totalPlusOneMesh,totalMinusOneMesh] = Amplitudes(t);
            set(p1,'CData', totalProductMesh);
            set(p2,'CData', totalPlusOneMesh);
            set(p3,'CData', totalMinusOneMesh);
            drawnow;
        end
    end

    function [xMesh,zMesh,totalProductMesh,totalPlusOneMesh,totalMinusOneMesh] = Amplitudes(t)
        % This would be much nicer done with a piecewise function...
        Xmain = 0:0.1:10;
        Zmain = 0:0.01:1;
        Xac = -5:0.1:-0.1;
        Zop = -0.5:0.01:-0.01;
        
        [xMeshMain,zMeshMain] = meshgrid(Xmain,Zmain);
        [xMeshOptic,zMeshOptic] = meshgrid(Xmain,Zop);
        [xMeshAcoustic,zMeshAcoustic] = meshgrid(Xac,Zmain);
        productMesh = Product(xMeshMain,zMeshMain,t);
        plusOneMesh = PlusOneAmp(xMeshMain,zMeshMain,t);
        minusOneMesh = MinusOneAmp(xMeshMain,zMeshMain,t);
        acousticMesh = AcousticAmp(xMeshAcoustic,zMeshAcoustic,t);
        opticMesh = OpticAmp(xMeshOptic,zMeshOptic,t);
        
        [xMesh,zMesh] = meshgrid([Xac Xmain],[Zop Zmain]);
        totalProductMesh = [ Zop'*Xac*0, opticMesh; acousticMesh, productMesh];
        totalPlusOneMesh = [ Zop'*Xac*0, opticMesh; acousticMesh, plusOneMesh];
        totalMinusOneMesh = [ Zop'*Xac*0, opticMesh; acousticMesh, minusOneMesh];
    end

    function amp = PlusOneAmp(x,z,t)
        amp = cos( (w_op + w_ac).*t - k_ac.*x - k_op.*z );
    end
    function amp = MinusOneAmp(x,z,t)
        amp = cos( (w_op - w_ac).*t + k_ac.*x - k_op.*z );
    end
    function amp = Product(x,z,t)
        amp = OpticAmp(x,z,t).*AcousticAmp(x,z,t);
    end
    function amp = OpticAmp(x,z,t)
        amp = cos( w_op.*t - k_op.*z );
    end
    function amp = AcousticAmp(x,z,t)
        amp = cos( w_ac.*t - k_ac.*x );
    end
end

