function [ opt ] = optimise_aod_orientation()

microSecs = -4:2:4;
xMilsArr = -20:10:20;
yMilsArr = xMilsArr;
[xMils, yMils] = meshgrid(xMilsArr, yMilsArr);
xMils = xMils(:)';
yMils = yMils(:)';
focalLength = 5;
xDeflectMils = 0;
yDeflectMils = xDeflectMils;
[xDeflectMils, yDeflectMils] = meshgrid(xDeflectMils, yDeflectMils);
xDeflectMils = xDeflectMils(:)';
yDeflectMils = yDeflectMils(:)';
pairDeflectionRatio = 0.4;
baseFreq = 30e6;
opt = Opt4range();

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

    function [opt] = Simple4()
        tic
        numberOfAods = 4;
        theta = cell(1,numberOfAods);
        phi = cell(1,numberOfAods);
        theta{1} = 0.023;
        phi{1} = -pi/2;
        theta{2} = 0.023;
        phi{2} = -pi/2;
        theta{3} = 0.025;
        phi{3} = -1.6176;
        theta{4} = 0;
        phi{4} = 1.6;
        
        [theta,phi] = GetAnglePermutations(numberOfAods, [theta,phi]);
        val = aol_performance(microSecs,xMils,yMils, theta, phi, xDeflectMils, yDeflectMils, pairDeflectionRatio, baseFreq, true,4 );
        PlotFovEfficiency(theta,phi)
        opt = [theta,phi];
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
        t = [0.02 0.02 0.04 0];
        p = [-pi/2   -pi/2   -2.7   -pi/2];
        for n = 1:4
            opt = Aod4(n, t, p);
            t(n) = opt(1);
            p(n) = opt(2);
        end
        eff = aol_performance(microSecs,xMils,yMils, t, p, xDeflectMils, yDeflectMils, pairDeflectionRatio, baseFreq, true,4 );
        PlotFovEfficiency(t,p);
        opt = [eff,t,p];
        
        function [opt] = Aod4(n,tTest,pTest)
            tic
            v = [tTest(n) pTest(n)];
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