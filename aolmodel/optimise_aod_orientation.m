function [ opt ] = optimise_aod_orientation()

microSecs = 4;
xyMm = GeneratePositionGrid(2);
xyDeflectMm = GeneratePositionGrid(5);
pairDeflectionRatio = 0;
baseFreq = 40e6;
opt = Simple4(0);

    function [opt] = Simple4(scanSpeed)
        tic
        t = [0.0389    0.0530    0.0029    0.0015 ];
        p = [-1.5708   -2.3157   -2.7    -1.5708]; % 30MHz
        aolPerturbations = aol_perturbations(t,p);
        t = [0.0399    0.0639   -0.0175   -0.0121];
        p = [ -1.5708  -2.4666   -2.7001   -1.5708]; % 40MHz
        %aolPerturbations.AddPerturbation(t,p);
        t = [0 0 0 0];
        p = [0 0 0 0];
        %aolPerturbations.AddPerturbation(t,p);
        
        val = aol_efficiency(microSecs,xyMm, aolPerturbations, xyDeflectMm, pairDeflectionRatio, baseFreq, scanSpeed, 4, true );
        opt = val;
        toc
    end

    function opt = Opt4range()
        opt = zeros(4,9);
        for m = 0:size(opt,1)
            pairDeflectionRatio = 1/3*m;
            opt(m+1,:) = Opt4();
        end
    end

    function [opt] = Opt4()
        t = [0.0389    0.0530    0.0029    0.0015]; %30MHz
        p = [-1.5708   -2.3157   -2.7    -1.5708]; % 30MHz
        t = [ 0.0399    0.0638   -0.0399   -0.0095]; %40MHz
        p = [ -1.5708  -2.4666    -0.2385   -1.5708]; % 40MHz
        for n = 1:4
            opt = Aod4(n, t, p);
            t(n) = opt(1);
            p(n) = opt(2);
        end
        eff = aol_efficiency(microSecs,xMm,yMm, t, p, xDeflectMm, yDeflectMm, pairDeflectionRatio, baseFreq, true,4,0 );
        opt = [eff,t,p];
        
        function [opt] = Aod4(n,tTest,pTest)
            tic
            v = [tTest(n) pTest(n)];
            opt = fminunc(@MinFun,v);
            toc
            function val = MinFun(v)
                tTest(n) = v(1);
                pTest(n) = v(2);
                val = -aol_efficiency(microSecs,xMm,yMm, tTest, pTest, xDeflectMm, yDeflectMm, pairDeflectionRatio, baseFreq, false, n,0 );
                val = harmmean(val,2);  % average over all deflections: want low efficiency in fov to have big effect
            end
        end
    end

    function PlotFovEfficiency(theta,phi)
        fovMm = focalLength * 100;
        divs = 10;
        [xDeflectMmLocal,yDeflectMmLocal] = meshgrid(linspace(-fovMm/2,fovMm/2,divs));
        xDeflectMmLocal = xDeflectMmLocal(:)';
        yDeflectMmLocal = yDeflectMmLocal(:)';
        [ deflectionEff ] =  aol_efficiency( microSecs, xMm, yMm, theta, phi, xDeflectMmLocal, yDeflectMmLocal, pairDeflectionRatio, baseFreq, false, size(theta,2));
        figure();
        contour(reshape(xDeflectMmLocal/focalLength,divs,divs),reshape(yDeflectMmLocal/focalLength,divs,divs),reshape(deflectionEff,divs,divs));
        colorbar;
        xlabel('x millirad')
        ylabel('y millirad')
    end
end

function xyMm = GeneratePositionGrid(gridSpacing) 
xMmArr = gridSpacing;
yMmArr = xMmArr;
[xMm, yMm] = meshgrid(xMmArr, yMmArr);
xMm = xMm(:)';
yMm = yMm(:)';
xyMm = [xMm;yMm];
end