function [ opt ] = optimise_aod_orientation()
tic

microSecs = -10:10:10;
xyMm = GeneratePositionGrid(-5:5:5);
xyDeflectMm = GeneratePositionGrid(10);%[0,10,20;10,10,10];% 
pairDeflectionRatio = 0;
baseFreq = 40e6;
scanSpeed = [000;000];
focalLength = 1;
aolPerturbations = SpecifyPerturbations(-1);
transducerWidths = [3.6 3.6 1.8 1.8] * 1e-3;
opWavelenVac = 800e-9;
   
opt = Simple4();

%opt = OptimiseAngles4();

%PlotFovEfficiencyPointing(pairDeflectionRatio,focalLength);
%PlotFovEfficiencyScanning(1,1,[200,600,1400,2000,3000,4000,5000; 0 0 0 0 0 0 0]);

    function [eff] = Simple4()
        driveParams =  aol_drive_params.MakeDriveParams(xyDeflectMm, pairDeflectionRatio, scanSpeed, baseFreq, focalLength, opWavelenVac);
        eff = plot_and_analyse_aol(microSecs,xyMm, aolPerturbations, driveParams, 4, transducerWidths, true );
    end

    function [opt] = OptimiseAngles4()
        t = [0.04    0.04   0   0];
        p = [ -pi/2  -1.5708   -1.5708   -1.5708];
        driveParamsForOptimising =  aol_drive_params.MakeDriveParams([0;0], 1, 0, baseFreq, 5, opWavelenVac);
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
        driveParams =  aol_drive_params.MakeDriveParams(xyDeflectMmLocal, pairDeflectionRatioLocal, 0, baseFreq, focalLengthLocal, opWavelenVac);
        deflectionEff = plot_and_analyse_aol(microSecs,xyMm, SpecifyPerturbations(-1), driveParams, 4, transducerWidths, true );
        figure();
        contour(xDeflectMmLocal/focalLengthLocal,yDeflectMmLocal/focalLengthLocal,reshape(deflectionEff,divs,divs));
        colorbar;
        xlabel('x millirad')
        ylabel('y millirad')
    end

    function PlotFovEfficiencyScanning(pairDeflectionRatioLocal, focalLengthLocal, scanSpeedLocal)
        scanSpeedLocal = -scanSpeedLocal;
        fovAngle = 36e-3;
        divs = 11;
        angles = linspace(-fovAngle,fovAngle,divs);
        
        xMms = focalLengthLocal * 1000 * angles;
        deflectionEff = zeros(divs,length(pairDeflectionRatioLocal)*length(baseFreq),length(scanSpeedLocal));
        for speedNum = 1:length(scanSpeedLocal)
            microSecsLocal = xMms / scanSpeedLocal(1,speedNum) * 1000;
            driveParams =  aol_drive_params.MakeDriveParams(GeneratePositionGrid(0), pairDeflectionRatioLocal, scanSpeedLocal(:,speedNum), baseFreq, focalLengthLocal, opWavelenVac);
            deflectionEff(:,:,speedNum) = plot_and_analyse_aol(microSecsLocal,xyMm, SpecifyPerturbations(-1), driveParams, 4, transducerWidths, false );
        end
        figure()
        plot(angles*1000,reshape(deflectionEff,divs,numel(deflectionEff)/divs))
        legend(cellstr(num2str(scanSpeedLocal(1,:)'/focalLengthLocal, '%-d')))
        xlabel('angle / mrad')
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