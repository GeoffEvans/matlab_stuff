function amplitude_modulation()
    transducer_width = 3.2e-3;

    figure; hold on;
    xlabel('Deviation from Bragg angle (degrees)')
    ylabel('Diffraction efficiency')
    rads = linspace(0, 23.52e-3);
    angular_deviation = rads/2/pi*360*2.226;
    powers = linspace(12, 0);
    max_effs = zeros(size(rads));
    max_effs_low = zeros(size(rads));
    for p = powers
        amp = amplitude(p, transducer_width);
        effs = get_eff(rads, amp, transducer_width);
        plot(angular_deviation, effs, 'Color', [(p > 2) && (p <= 4) || (p > 8) && (p <= 10), (p > 4) && (p <= 8), (p < 2) || (p > 6) && (p <= 8) || (p > 2) && (p <= 4)]);
        max_effs = max(effs, max_effs);
        if p < 4
            max_effs_low = max(effs, max_effs_low);
        end
    end
    set(gca,'TickDir','out')
    figure()
    hold on
    plot1 = plot(angular_deviation, max_effs, 'r');
    plot2 = plot(angular_deviation, get_eff(rads, amplitude(10.5, transducer_width), transducer_width), 'r-');
    plot3 = plot(angular_deviation, max_effs_low, 'b');
    plot4 = plot(angular_deviation, get_eff(rads, amplitude(1.25, transducer_width), transducer_width), 'b-');
    set(plot1,'Color',[0.6 0 0], 'linewidth', 2);
    set(plot2,'Color',[1 0.4 0.4], 'linewidth', 2);
    set(plot3,'Color',[0 0 1], 'linewidth', 2);
    set(plot4,'Color',[0.3 0.75 0.93], 'linewidth', 2);
    set(gca,'TickDir','out')
end

    
function eff = get_eff(angles, amp, transducer_width)
    
    wavelength = 800e-9;
    k = 2*pi / wavelength;
    wavevector_mismatches_mag = 160e3 * angles;
    
    n_in = 2.2264;
    n_out = 2.2262;
    p = -0.12;

    delta_n0 = -0.5 * power(n_in, 2.) * n_out * p * amp;
    delta_n1 = -0.5 * power(n_out, 2.) * n_in * p * amp;
    v0 = - k * delta_n0 * transducer_width;
    v1 = - k * delta_n1 * transducer_width;

    zeta = -0.5 * wavevector_mismatches_mag * transducer_width;
    sigma = sqrt(power(zeta, 2.) + v0.*v1./4);
    
    sinc = sin(sigma) ./ sigma;
    eff = v0.*v1./4 .* power(sinc, 2.);
end
function amp = amplitude(power, transducer_width)
    teo2_density = 5990;
    transducer_height = 15e-3;
    vel = 613;
    numerator = 2 * power;
    denominator = teo2_density * vel^3 * transducer_width * transducer_height;
    amp = sqrt(numerator ./ denominator);
end
