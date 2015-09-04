function psf_4()
    time = 0e-6;
    w2 = 0
    w4 = 0
    ws = 0
    x1 = [0, w2, 0, w4, 0];
    y1 = [0, w2, 0, w4, 0];    
    x2 = [0, w2, 0, w4, 0];
    y2 = [0, w2, 0, w4, 0];    
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    a = aol_fft();
    a.z_list = a.z_list + 10e-6;
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    a.xy_xz_plot_psf(propagated_wave, x, y, 0.5);
    res = a.get_psf_dimensions(propagated_wave, x, y, a.z_list);
end

function ideal_correction()
tic
    waves = 12;
    scale = 3;
    list = -waves:waves;
    results = zeros(11, 3);
    for n = list
        results(n+waves+1, :) = run(9, 0, n/scale);
    end
    figure()
    hold on
    plot(list/scale, results(:,1), 'r')
    plot(list/scale, results(:,2), 'b')
    plot(list/scale, results(:,3), 'g')
    hold off
    
    xlabel('waves of spherical')
    ylabel('(red) x FWHM (blue) z FWHM (green) max intensity squared')
toc
end

function aod_correction()
tic
    waves = 12;
    scale = 4;
    list = -waves:waves;
    results = zeros(11, 3);
    for n = list
        results(n+waves+1, :) = run(9, n/scale, 0);
    end
    figure()
    hold on
    plot(list/scale, results(:,1), 'r')
    plot(list/scale, results(:,2), 'b')
    plot(list/scale, results(:,3), 'g')
    hold off
    
    xlabel('waves of d')
    ylabel('(red) x FWHM (blue)z FWHM (green) max intensity squared')
toc
end

function focus_effect()
    waves = 10;
    results = zeros(11, 3);
    for n = -waves:waves
        results(n+waves+1, :) = run(n, 0, 0);
    end
    figure()
    hold on
    plot(-waves:waves, results(:,1), 'r')
    plot(-waves:waves, results(:,2), 'b')
    plot(-waves:waves, results(:,3), 'g')
    hold off
    
    xlabel('waves of focus')
    ylabel('(red) x FWHM (blue)z FWHM (green) max intensity squared')
end

function res = run(w2, w4, ws)
    time = 0e-6;
    x1 = [0, w2, 0, w4, 0];
    y1 = [0, w2, 0, w4, 0];    
    x2 = [0, w2, 0, w4, 0];
    y2 = [0, w2, 0, w4, 0];    
    
    waves = struct();
    waves.focus = 0;
    waves.spherical = ws;
    waves.ac_angle = 0; % walk off angle in degrees
    waves.aods = {x1, y1, x2, y2};
    
    a = aol_fft();
    [propagated_wave, x, y] = a.calculate_psf_through_aol(waves, time);
    %a.xy_xz_plot_psf(propagated_wave, x, y, 0.5);
    res = get_psf_dimensions(propagated_wave, x, y, a.z_list);
end
