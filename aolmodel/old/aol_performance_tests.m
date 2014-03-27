function aol_performance_tests()

microSecs = -10:10:10;
xMilsArr = -20:20:20;
yMilsArr = xMilsArr;
[xMils, yMils] = meshgrid(xMilsArr, yMilsArr);
xMils = xMils(:)';
yMils = yMils(:)';
xDeflectMils = -0;
yDeflectMils = xDeflectMils*2;
[xDeflectMils, yDeflectMils] = meshgrid(xDeflectMils, yDeflectMils);
xDeflectMils = xDeflectMils(:)';
yDeflectMils = yDeflectMils(:)';
pairDeflectionRatio = 2;

pass = TestRaysPassStraightThroughAtZeroFreq();
fprintf([num2str(pass) '\t TestRaysPassStraightThroughAtZeroFreq \n' ])

pass = TestRaysGoToCorrectPoint();
fprintf([num2str(pass) '\t TestRaysGoToCorrectPoint \n' ])

pass = TestRaysFocus();
fprintf([num2str(pass) '\t TestRaysFocus \n' ])

    function pass = TestRaysPassStraightThroughAtZeroFreq()
        [ ~,x,y ] =  aol_performance( microSecs, 0,0, [0 0 0 0], [0 0 0 0], 0, 0, pairDeflectionRatio, 0, false, 1 );
        xPass = max(abs(x(:))) < 1e-14;
        yPass = max(abs(y(:))) < 1e-14;
        pass = xPass && yPass;
    end
    function pass = TestRaysGoToCorrectPoint()
        theta = [0 0 0 0;0.06 0.07 0.08 0.04];
        phi = [0 0 0 0;3 2 1 2];
        baseFreq = 35e6;
        [ ~,x,y ] =  aol_performance( microSecs, xMils,yMils, theta, phi, 10, -4, pairDeflectionRatio, baseFreq, false, 1 );
        xOff = x(end-1,:) - 10e-3;
        yOff = y(end-1,:) + 4e-3;
        xPass = max(abs(xOff)/x(end-1,:)) < 1e-1;
        yPass = max(abs(yOff)/y(end-1,:)) < 1e-1;
        pass = xPass && yPass;
    end
    function pass = TestRaysFocus()
        theta = [0 0 0 0;0.06 -0.06 0.06 0.06];
        phi = [0 0 0 0;3 2 1 2];
        baseFreq = 35e6;
        numOfTimes = length(microSecs);
        numOfPositions = length(xMils);
        numOfDeflections = 1;
        numOfPerturbations = size(theta,1);
        [ ~,x,y ] =  aol_performance( microSecs, xMils,yMils, theta, phi, 5, 17, pairDeflectionRatio, baseFreq, false, 1 );
        xFocus = x(end-2,:);
        yFocus = y(end-2,:);
        sigmaX = std(reshape(xFocus,numOfTimes*numOfPositions,numOfDeflections*numOfPerturbations),1);
        sigmaY = std(reshape(yFocus,numOfTimes*numOfPositions,numOfDeflections*numOfPerturbations),1);
        xPass = sigmaX < 1e-4;
        yPass = sigmaY < 1e-4;
        pass = min(xPass & yPass);
    end
end