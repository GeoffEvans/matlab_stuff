function [ opt ] = optimise_aod_orientation_()

microSecs = -4:4:4;
xMilsArr = -20:20:20;
yMilsArr = xMilsArr;
[xMils, yMils] = meshgrid(xMilsArr, yMilsArr);
xMils = xMils(:)';
yMils = yMils(:)';
focalLength = 2;
xDeflectMils = 0;
yDeflectMils = xDeflectMils;
[xDeflectMils, yDeflectMils] = meshgrid(xDeflectMils, yDeflectMils);
xDeflectMils = xDeflectMils(:)';
yDeflectMils = yDeflectMils(:)';
pairDeflectionRatio = 0.4;
baseFreq = 40e6;
opt = Simple4(0);

    function [opt] = Simple2()
        tic
        numberOfAods = 2;
        theta = cell(1,numberOfAods);
        phi = cell(1,numberOfAods);
        theta{1} = 0.023;
        phi{1} = -pi/2;
        theta{2} = 0 / 180 * pi;
        phi{2} = 90 / 180 * pi;
        
        [theta,phi] = GetAnglePermutations(numberOfAods, [theta,phi]);
        val = -aol_performance(microSecs,xMils,yMils, theta, phi, xDeflectMils, yDeflectMils, pairDeflectionRatio, baseFreq, true, 2 );
        PlotFovEfficiency(theta,phi)
        opt = [theta,phi];
        toc
    end

    function [opt] = Opt2()
        t = [0.0392    0.0942   ];
        p = [4.7124    0.0051   ];
        for n = 1:2
            opt = Aod2(n, t, p);
            t(n) = opt(1);
            p(n) = opt(2);
        end
        PlotFovEfficiency(theta,phi)
        opt = [t,p];
        
        function [opt] = Aod2(n,tTest,pTest)
            tic
            v = [tTest(n) pTest(n)];
            lb = [0 -pi];
            ub = [0.1 pi];
            display('opt')
            opt = fminunc(@MinFun,v);
            toc
            function val = MinFun(v)
                tTest(n) = v(1);
                pTest(n) = v(2);
                val = -aol_performance(microSecs,xMils,yMils, tTest, pTest, xDeflectMils, yDeflectMils, pairDeflectionRatio, baseFreq, false, n );
                val = harmmean(val,2);  % average over all deflections: want low efficiency in fov to have big effect
            end
        end
    end

    function [opt] = Simple4(scanSpeed)
        tic
        t = [0.0389    0.0530    0.0029    0.0015]; %30MHz
        p = [-1.5708   -2.3157   -2.7    -1.5708]; % 30MHz
        t = [0.0399    0.0639   -0.0175   -0.0121]; %40MHz
        p = [ -1.5708  -2.4666   -2.7001   -1.5708]; % 40MHz
        val = aol_performance(microSecs,xMils,yMils, t, p, xDeflectMils, yDeflectMils, pairDeflectionRatio, baseFreq, true,4, -scanSpeed );
        %PlotFovEfficiency(t,p)
        opt = val;
        toc
    end

    function d = CompareScans()
        d = zeros(7,13);
        scanSpeed = 72/48*1000;
        microSecs = -24:4:24;
        d(1,:) = Simple4(scanSpeed);
        
        scanSpeed = 72/60*1000;
        microSecs = -30:5:30;
        d(2,:) = Simple4(scanSpeed);

        scanSpeed = 72/90*1000;
        microSecs = -45:7.5:45;
        d(3,:) = Simple4(scanSpeed);
        
        scanSpeed = 72/120*1000;
        microSecs = -60:10:60;
        d(4,:) = Simple4(scanSpeed);
        
        scanSpeed = 72/180*1000;
        microSecs = -90:15:90;
        d(5,:) = Simple4(scanSpeed);
        
        scanSpeed = 72/360*1000;
        microSecs = -180:30:180;
        d(6,:) = Simple4(scanSpeed);
        
        scanSpeed = 0;
        microSecs = -4:1:4;
        for xDeflectMils = -36:6:36
            d(7,-xDeflectMils/6+7) = mean(Simple4(0));
        end
        
        figure()
        plot(-(-18:3:18)'*ones(1,7),d');
        legend('1500','1200','800','600','400','200','pointing');
        xlabel('deflection / mrad')
        ylabel('efficiency')
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
        eff = aol_performance(microSecs,xMils,yMils, t, p, xDeflectMils, yDeflectMils, pairDeflectionRatio, baseFreq, true,4,0 );
        %PlotFovEfficiency(t,p);
        opt = [eff,t,p];
        
        function [opt] = Aod4(n,tTest,pTest)
            tic
            v = [tTest(n) pTest(n)];
            opt = fminunc(@MinFun,v);
            toc
            function val = MinFun(v)
                tTest(n) = v(1);
                pTest(n) = v(2);
                val = -aol_performance(microSecs,xMils,yMils, tTest, pTest, xDeflectMils, yDeflectMils, pairDeflectionRatio, baseFreq, false, n,0 );
                val = harmmean(val,2);  % average over all deflections: want low efficiency in fov to have big effect
            end
        end
    end

    function PlotFovEfficiency(theta,phi)
        fovMils = focalLength * 100;
        divs = 10;
        [xDeflectMilsLocal,yDeflectMilsLocal] = meshgrid(linspace(-fovMils/2,fovMils/2,divs));
        xDeflectMilsLocal = xDeflectMilsLocal(:)';
        yDeflectMilsLocal = yDeflectMilsLocal(:)';
        [ deflectionEff ] =  aol_performance( microSecs, xMils, yMils, theta, phi, xDeflectMilsLocal, yDeflectMilsLocal, pairDeflectionRatio, baseFreq, false, size(theta,2));
        figure();
        contour(reshape(xDeflectMilsLocal/focalLength,divs,divs),reshape(yDeflectMilsLocal/focalLength,divs,divs),reshape(deflectionEff,divs,divs));
        colorbar;
        xlabel('x millirad')
        ylabel('y millirad')
    end
end

function [theta,phi] = GetAnglePermutations(numberOfAods,sets)
cartProd = CartesianProduct(sets);
theta = cartProd(:,1:numberOfAods);
phi = cartProd(:,numberOfAods+1:numberOfAods*2);

    function result = CartesianProduct(sets)
        c = cell(1, numel(sets));
        [c{:}] = ndgrid( sets{:} );
        result = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
    end
end