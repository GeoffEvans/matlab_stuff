function [ opt ] = optimise_aod_orientation()

microSecs = -4:4:4;
xMilsArr = -20:40:20;
yMilsArr = xMilsArr;
[xMils, yMils] = meshgrid(xMilsArr, yMilsArr);
xMils = xMils(:)';
yMils = yMils(:)';
xDeflectMils = -30:60:30;
yDeflectMils = xDeflectMils*2;
[xDeflectMils, yDeflectMils] = meshgrid(xDeflectMils, yDeflectMils);
xDeflectMils = xDeflectMils(:)';
yDeflectMils = yDeflectMils(:)';
pairDeflectionRatio = 2;
opt = Brute();

    function [opt] = Brute()  
        tic
        numberOfAods = 4;
        theta = cell(1,numberOfAods);
        phi = cell(1,numberOfAods);
        theta{1} = 0 / 180 * pi;
        phi{1} = 270 / 180 * pi;
        theta{2} = 0 / 180 * pi;
        phi{2} = 90 / 180 * pi;
        theta{3} = 0 / 180 * pi;
        phi{3} = 270 / 180 * pi;
        theta{4} = 0 / 180 * pi;
        phi{4} = 90 / 180 * pi;
        sets = [theta,phi];
        
        [theta,phi] = GetAnglePermutations(numberOfAods, sets);
        [ prodEffOpt ] =  aol_performance( microSecs, xMils,yMils, theta, phi, xDeflectMils, yDeflectMils, pairDeflectionRatio, 30e6, true, 1 );
        opt = prodEffOpt;
        toc
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
    end

    function [opt] = Aod2()
        tic
        MinFun = @(v) -aol_performance(microSecs,xMils,yMils, [2.1990/180*pi, v(1)], [270/180*pi, v(2)], false, 2);
        v = [0, 90] / 180 * pi;
        v = fminsearch(MinFun,v);
        optEff = aol_performance(microSecs,xMils,yMils, [2.1990/180*pi, v(1)], [270/180*pi, v(1)], true, 2);
        opt = [optEff v*180/pi];
        toc
    end
    function [opt] = Opt2()
        t = [2.1990, 0.0446] / 180 * pi; 
        p = [270, 89.9881] / 180 * pi;  
        opt = aol_performance(microSecs,xMils,yMils, t, p, false, 2);
    end

    function [opt] = Opt4()
        t = [0.0385    0   -0.0021    0];
        p = [4.7124    0    2.8315    0];
        for n = 1:4
            opt = Aod4(n, t, p);
            t(n) = opt(1);
            p(n) = opt(2);
        end        
        opt = aol_performance(microSecs,xMils,yMils, t, p, xDeflectMils, yDeflectMils, pairDeflectionRatio, 30e6, true, 4);
        opt = [opt t p];
        
        function [opt] = Aod4(n,t,p)
            tic
            v = [t(n) p(n)];
            lb = [-1 -2*pi];
            ub = [1 2*pi];
            opt = simulannealbnd(@MinFun,v,lb,ub);
            toc
            
            function val = MinFun(v) 
                t(n) = v(1);
                p(n) = v(2);
                val = -aol_performance(microSecs,xMils,yMils, t, p, xDeflectMils, yDeflectMils, pairDeflectionRatio, 30e6, false, n );
            end
            %  0.8859    0.0407    0.0408   -0.0087   -0.0292    4.7124     3.1482    5.9256    3.1074
            %  0.8859    0.0407    0.0407   -0.0028   -0.0302    4.7124     3.1482    4.7124    3.1063      
            %  0.8859    0.0407    0.0396    0.0101    0.0133    4.7124     3.3055    2.8308    2.3055
        end
    end
end