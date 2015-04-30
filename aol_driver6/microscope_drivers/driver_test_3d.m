function driver_test_3d()
    wavelength = 920e-9;
    optimalFrequency = ones(1,6) * 39 * 1e6;
    acceptanceAngle = 4; % mrad
    zoomFactor = 1;
    xyNumOfElems = 3;
    dwellTime = 10e-6 * 1;
    imageCentreNormalised = [[0;0;1]];% [-3.5;0;1], [-2;1;1.2]]; 
    scanDispNormalised = [[0;0;0]];% [-0.2;0;0], [-0.2;-1;0]];
    aodAperture = 15e-3;
    V = 612.9;
    aodMode = -1;
    imagingMode = 'structural';

    tic
    [baseFreqs, linearChirps, rampTime, isMiniscan] = microscope_driver6(...
        imagingMode,...    
        aodMode,...
        xyNumOfElems,...
        acceptanceAngle,... % mrad
        dwellTime,... 
        zoomFactor,...
        optimalFrequency,...
        imageCentreNormalised,...
        scanDispNormalised,...
        aodAperture,...
        V,...
        wavelength); 
    toc
    
    %baseFreqs = baseFreqs*0;
    %linearChirps = Aod6cyclicNewScan(0.8)';
    numOfRays = 11;
    numOfDrives = size(baseFreqs,1);
    spacings = [ones(1,5)*5e-2, 0.9375];
    acDirections = 0.5 * [2,1,-1,-2,-1,1;0,sqrt(3),sqrt(3),0,-sqrt(3),-sqrt(3)];

    figure();
    hold on
    for n = -2:2
        t = n * rampTime / 2;
        if isMiniscan
            t = ones(numOfRays,1)*t;
        end

        kx = zeros(numOfRays,numOfDrives);
        ky = zeros(numOfRays,numOfDrives);
        x = ( cumsum(ones(numOfRays,numOfDrives), 1) - (numOfRays+1)/2 ) * 1e-2;
        y = x;
        z = 0*x - sum(spacings(1:end-1));
        for m = 1:6
            spacing = spacings(m);
            acDir = acDirections(:,m);
            centre = [0,0];
            baseFreq = repmat(baseFreqs(:,m)',numOfRays,1);
            linearChirp = repmat(linearChirps(:,m)',numOfRays,1);

            timeFromTransducerToBaseRay = ((x - centre(1))*acDir(1) + (y - centre(2))*acDir(2) + aodAperture/2) / V;
            T = (t - timeFromTransducerToBaseRay); % takes longer to get there so originated earlier
            aodFreq = baseFreq + linearChirp .* T;
            kx = kx - (wavelength / V) * acDir(1) * aodFreq;
            ky = ky - (wavelength / V) * acDir(2) * aodFreq;
            x_prev = x;
            x = x - kx .* spacing;
            y_prev = y;
            y = y - ky .* spacing;
            z_prev = z;
            z = z + spacing;
            for p = 1:numOfRays
                for q = 1:numOfDrives
                    plot([y_prev(p,q), y(p,q)], [z_prev(p,q), z(p,q)], 'color', (n + 2)/5 * [1, -1, 0] + [0, 1, 0]);
                end
            end
        end   
    end
    hold off     
end

function chirp = Aod6cyclicNewScan(focalLength) % z = 1/3f for second solution
    f = focalLength;
    chirpFactor = [... 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              (5*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/((2400*f^2 + 440*f + 43/3)*(4800*f^2 + 930*f + 127/3)) - (10*(16800*f^2 + 1700*f + 95/3))/((2400*f^2 + 440*f + 43/3)*(4800*f^2 + 930*f + 127/3));...
                  ((149140*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(9*(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9)) - 91200*f - 896000*f^2 + (28912000*f^2*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (284160000*f^3*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (921600000*f^4*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (3525200*f*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(3*(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9)) - 16360/9)/(46080000*f^4 + 16128000*f^3 + 1996000*f^2 + 105060*f + 1984);...
 -((149140*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(9*(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9)) - 408800*f - 5216000*f^2 - 19200000*f^3 + (28912000*f^2*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (284160000*f^3*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (921600000*f^4*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (3525200*f*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(3*(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9)) - 85520/9)/(11520000*f^4 + 3840000*f^3 + 481600*f^2 + 26880*f + 1700/3);...
                ((149140*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(9*(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9)) - 117600*f - 928000*f^2 + (28912000*f^2*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (284160000*f^3*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (921600000*f^4*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9) + (3525200*f*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(3*(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9)) - 31720/9)/(46080000*f^4 + 14688000*f^3 + 1676000*f^2 + 80560*f + 4064/3);...
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (20*(65000*f - (530841600000000*f^8 + 362741760000000*f^7 + 106131456000000*f^6 + 17353728000000*f^5 + 1733150720000*f^4 + (324409088000*f^3)/3 + (36986636800*f^2)/9 + (781360480*f)/9 + 63134449/81)^(1/2) + 1166400*f^2 + 8640000*f^3 + 23040000*f^4 + 3673/3))/(46080000*f^4 + 14208000*f^3 + 1445600*f^2 + (176260*f)/3 + 7457/9);...
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         1/(3*f)];
    chirp = 612.9^2 / 920e-9 * chirpFactor;
end