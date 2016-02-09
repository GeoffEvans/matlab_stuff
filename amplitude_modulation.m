function amplitude_modulation()
    transducer_width = 3.2e-3;

    figure; hold on;
    xlabel('incidence angle')
    ylabel('efficiency')
    rads = linspace(0, 20e-3);
    powers = linspace(8, 0);
    for p = powers
        amp = amplitude(p, transducer_width);
        effs = get_eff(rads, amp, transducer_width);
        plot(rads/2/pi*360*2.226, effs, 'Color', [(p > 2) && (p < 4), (p > 4) && (p < 6), (p < 2)]);
    end
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
