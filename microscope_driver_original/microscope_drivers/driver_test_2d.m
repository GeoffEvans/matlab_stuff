function driver_test_2d()

wavelength = 920 * 1e-9;
optimalFrequency = [1,1,1,1] * 39 * 1e6;
acceptanceAngle = 18; % mrad
zoomFactor = [1];
xyNumOfElems = 5;
dwellTime = 10e-6 * 1;
imageCentreNormalised = [[0;0;0]]% [-3.5;0;1], [-2;1;1.2]]; 
scanDispNormalised = [[0;0;0]]% [-0.2;0;0], [-0.2;-1;0]];
pairDeflectionRatio = -1;
aodAperture = 15e-3;
V = 612.8;
aodMode = -1;
imagingMode = 'structural';

tic
[baseFreq, linearChirp, rampTime, isMiniscan] = microscope_driver(...
    imagingMode,...    
    aodMode,...
    xyNumOfElems,...
    acceptanceAngle,... % mrad
    dwellTime,... 
    zoomFactor,...
    optimalFrequency,...
    imageCentreNormalised,...
    scanDispNormalised,...
    pairDeflectionRatio,...
    aodAperture,...
    V,...
    wavelength); 
toc

V = teo2.find_v_ac_min(pi/2,pi/4);
aodAperture = 15e-3;
numOfRays = 11;
numOfDrives = size(baseFreq,2);

aod1 = struct(...
    'baseFreq',repmat(baseFreq(1,:),numOfRays,1),...
    'linearChirp',repmat(linearChirp(1,:),numOfRays,1),...
    'l',0.1,...
    'transducerPosition',0-aodAperture/2);
aod2 = struct(...
    'baseFreq',repmat(baseFreq(3,:),numOfRays,1),...
    'linearChirp',repmat(linearChirp(3,:),numOfRays,1),...
    'l',0.05+0.1,...
    'transducerPosition',-4.6e-3+aodAperture/2);
lens = struct(...
    'focalLength', 0.1,...
    'l', 0.20e0,...
    'centre',-0.00522122);

figure();
for n = -2:2
    t = n * rampTime / 2;
    if isMiniscan
        t = ones(numOfRays,1)*t;
    end
    k = zeros(numOfRays,numOfDrives);
    x1 = ( cumsum(ones(numOfRays,numOfDrives)) - (numOfRays+1)/2 ) * 1e-2;

    timeFromTransducerToBaseRay1 = (x1 - aod1.transducerPosition) / V;
    T = (t - timeFromTransducerToBaseRay1);
    aodFreq1 = aod1.baseFreq + aod1.linearChirp .* T;
    k = k + (wavelength / V) * (-1) * aodFreq1;
    x2 = x1 + k .* aod1.l;

    timeFromTransducerToBaseRay2 = - (x2 - aod2.transducerPosition) / V;
    T = (t - timeFromTransducerToBaseRay2);
    aodFreq2 = aod2.baseFreq + aod2.linearChirp .* T;
    k = k + (wavelength / V) .* aodFreq2;
    x3 = x2 + k .* aod2.l;

    k = k - (x3 - lens.centre)/lens.focalLength;
    x4 = x3 + k .* lens.l;
    
    hold on;
    d = 0.1;
    z = cumsum([0,d,aod1.l,aod2.l,lens.l]);
    for p = 1:numOfRays
        for q = 1:numOfDrives
            plot(z,[x1(p,q),x1(p,q),x2(p,q),x3(p,q),x4(p,q)], 'color', (n + 2)/5 * [1, -1, 0] + [0, 1, 0]);
        end
    end
    hold off;
end
end