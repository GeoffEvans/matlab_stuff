function [ opt ] = optimise_aod_orientation()

microSecs = -4:1:4;
xMilsArr = 0;%-5:1:5;
yMilsArr = 0;%-15:5:5;
[xMils, yMils] = meshgrid(xMilsArr, yMilsArr);
xMils = xMils(:)';
yMils = yMils(:)';
opt = Brute();

    function [opt] = Brute()  
        tic
        numberOfAods = 4;
        theta = cell(1,numberOfAods);
        phi = cell(1,numberOfAods);
        theta{1} = 2.2 / 180 * pi;
        phi{1} = 270 / 180 * pi;
        theta{2} = 0 / 180 * pi;
        phi{2} = 90 / 180 * pi;
        theta{3} = 1.2 / 180 * pi;
        phi{3} = 270 / 180 * pi;
        theta{4} = 0 / 180 * pi;
        phi{4} = [90 270] / 180 * pi;
        sets = [theta,phi];
        
        [theta,phi] = GetAnglePermutations(numberOfAods, sets);
        [ prodEffOpt ] =  aol_performance( microSecs,xMils,yMils, theta, phi, true );
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

    function [opt] = Aod2Opt()
        tic
        MinFun = @(v) -aol_performance(microSecs,xMils,yMils, v(1:2), v(3:4), false);
        v = [2, 2, 270, 180] / 180 * pi;
        angles = fminsearch(MinFun,v);
        optEff = aol_performance(microSecs,xMils,yMils, angles(1:2), angles(3:4), true);
        opt = [optEff angles*180/pi];
        toc
    end

    function [opt] = Aod4Opt()
        tic
        MinFun = @(v) -aol_performance(microSecs,xMils,yMils, v(1:4), v(5:8), false);
        v = [2.2,2.3,0,0,270,190,90,0] / 180 * pi;
        options = optimset('MaxFunEvals', 2000);
        angles = fminsearch(MinFun,v,options);
        optEff = aol_performance(microSecs,xMils,yMils, angles(1:4), angles(5:8), true);
        opt = [optEff angles*180/pi];
        toc
    end
end