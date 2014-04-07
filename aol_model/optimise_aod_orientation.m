function [ opt ] = optimise_aod_orientation()
tic

microSecs = -4:4:4;
xyMm = GeneratePositionGrid(-5:5:5);
xyDeflectMm = GeneratePositionGrid(0);
pairDeflectionRatio = 1;
baseFreq = 40e6;
scanSpeed = 0;
focalLength = 2;
aolPerturbations = SpecifyPerturbations(-1);
transducerWidths = [3.6 3.6 1.8 1.8] * 1e-3;

%opt = Simple4();

%opt = OptimiseAngles4();

PlotFovEfficiencyPointing(pairDeflectionRatio,focalLength);
%PlotFovEfficiencyScanning(0,focalLength,[200,600,1400,2000,3000,4000,5000]);

    function [eff] = Simple4()
        driveParams = MakeDriveParams(xyDeflectMm, pairDeflectionRatio, scanSpeed, baseFreq, focalLength);
        eff = plot_and_analyse_aol(microSecs,xyMm, aolPerturbations, driveParams, 4, transducerWidths, true );
        if numel(eff) > 200
            display('more than 200 effs, returning max only')
            eff = max(eff(:));
        end
    end

    function [opt] = OptimiseAngles4()
        t = [0.04    0.04   0   0];
        p = [ -pi/2  -1.5708   -1.5708   -1.5708];
        driveParamsForOptimising = MakeDriveParams([0;0], 1, 0, baseFreq, 5);
        for nthAod = 1:4
            opt = FindOptimalPerturbation(nthAod, t, p);
            t(nthAod) = opt(1);
            p(nthAod) = opt(2);
        end
        eff = plot_and_analyse_aol(microSecs,xyMm, aol_perturbations(t,p), driveParamsForOptimising, 4, transducerWidths, true );
        opt = [eff,t,p];
        
        function [opt] = FindOptimalPerturbation(n,tTest,pTest)
            v = [tTest(n) pTest(n)];
            opt = simulannealbnd(@MinFun,v,[0 -pi],[0.1 pi])  % display each aod optimal
           
            function val = MinFun(v)
                tTest(n) = v(1);
                pTest(n) = v(2);
                val = -plot_and_analyse_aol(microSecs,xyMm, aol_perturbations(tTest,pTest), driveParamsForOptimising, n, transducerWidths, false );
            end
        end
    end

    function PlotFovEfficiencyPointing(pairDeflectionRatioLocal, focalLengthLocal)
        fovAngle = 36e-3;
        divs = 11;
        fovMm = focalLengthLocal * 1000 * fovAngle;
        [xDeflectMmLocal,yDeflectMmLocal] = meshgrid(linspace(-fovMm/2,fovMm/2,divs));
        xyDeflectMmLocal = [xDeflectMmLocal(:)'; yDeflectMmLocal(:)'];
        driveParams = MakeDriveParams(xyDeflectMmLocal, pairDeflectionRatioLocal, 0, baseFreq, focalLengthLocal);
        deflectionEff = plot_and_analyse_aol(microSecs,xyMm, SpecifyPerturbations(-1), driveParams, 4, transducerWidths, true );
        figure();
        contour(xDeflectMmLocal/focalLengthLocal,yDeflectMmLocal/focalLengthLocal,reshape(deflectionEff,divs,divs));
        colorbar;
        xlabel('x millirad')
        ylabel('y millirad')
    end

    function PlotFovEfficiencyScanning(pairDeflectionRatioLocal, focalLengthLocal, scanSpeedLocal)
        %scanSpeedLocal = -scanSpeedLocal;
        fovAngle = 36e-3;
        divs = 11;
        angles = linspace(-fovAngle,fovAngle,divs);
        
        xMms = focalLengthLocal * 1000 * angles;
        deflectionEff = zeros(divs,length(pairDeflectionRatioLocal)*length(baseFreq),length(scanSpeedLocal));
        for speedNum = 1:length(scanSpeedLocal)
            microSecsLocal = xMms / scanSpeedLocal(speedNum) * 1000;
            driveParams = MakeDriveParams(GeneratePositionGrid(0), pairDeflectionRatioLocal, scanSpeedLocal(speedNum), baseFreq, focalLengthLocal);
            deflectionEff(:,:,speedNum) = plot_and_analyse_aol(microSecsLocal,xyMm, SpecifyPerturbations(-1), driveParams, 4, transducerWidths, false );
        end
        figure()
        plot(angles*1000,reshape(deflectionEff,divs,numel(deflectionEff)/divs))
        legend(cellstr(num2str(scanSpeedLocal'/focalLengthLocal, '%-d')))
        xlabel('angle')
        ylabel('eff')
    end
toc
end

function xyMm = GeneratePositionGrid(gridSpacing)
    xMmArr = gridSpacing;
    yMmArr = xMmArr;
    [xMm, yMm] = meshgrid(xMmArr, yMmArr);
    xMm = xMm(:)';
    yMm = yMm(:)';
    xyMm = [xMm;yMm];
end

function driveParams = MakeDriveParams(xyDef,ratio,speed,optimalBaseFreq,focalLength)
    % returns arrays of horizontal form pure stretch on speed, pure repeat on ratio
    [xDefStretch,ratioStretched,speedStretched] = meshgrid(xyDef(1,:),ratio,speed);
    [yDefStretch,~,~] = meshgrid(xyDef(2,:),ratio,speed);
    xyDefStretched(1,:) = xDefStretch(:); % [ def1 x numOfRatios, def2 x numOfRatios, ...] x numOfSpeeds
    xyDefStretched(2,:) = yDefStretch(:);
    ratioStretched = ratioStretched(:)'; % repeat x numOfSpeedsByDefs
    speedStretched = speedStretched(:)'; % [ speed1 x numOfRatiosByDefs, speed2 x numOfRatiosByDefs, ... ]
    driveParams = aol_drive_params(focalLength, optimalBaseFreq, xyDefStretched, ratioStretched, speedStretched);
end

function aolPerturbations = SpecifyPerturbations(id)
    t = [0 0 0 0];
    p = [0 0 0 0];
    aolPerturbations = aol_perturbations(t,p);
    if ismember(1,id)
        t = [0.0399    0.0639   -0.0175   -0.0121];
        p = [ -1.5708  -2.4666   -2.7001   -1.5708]; % 40MHz
        aolPerturbations.AddPerturbation(t,p);
    end
    if ismember(2,id)
        t = [0.0389    0.0530    0.0029    0.0015 ];
        p = [-1.5708   -2.3157   -2.7    -1.5708]; % 30MHz
        aolPerturbations.AddPerturbation(t,p);
    end
    if ismember(-1,id)
        t = [ 0.0400    0.0644    0.0442    0.0117];
        p = [  -1.5708   -2.4695   2.8915    1.9074]; % 40MHz
        aolPerturbations = aol_perturbations(t,p);
    end
end