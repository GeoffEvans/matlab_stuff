function events = event_detection( data, duration, rise_time, decay_time, threshold, make_plot )

    if nargin == 4
        make_plot = false;
    end

    N = numel(data);
    step_time = duration/N;
    detection_criterion = zeros(1,N);
    times = 1:N;
    
    for time = times
        template = get_template(time, N, step_time, rise_time, decay_time);
        scale = get_scale(template, data, N);
        offset = (sum(data) - scale .* sum(template))/N;

        fitted_template = template .* scale + offset;
        sse = sum((data - fitted_template).^2);
        standard_error = sqrt(sse * (N-1));

        detection_criterion(time) = scale ./ standard_error;

    end
    
    events = times(get_events(threshold, detection_criterion));
    
    if make_plot
        hold on
        plot(times, detection_criterion)
        plot(events, events * 0, 'rx')
        hold off
    end
end

function events = get_events(threshold, criterion)
    thres = abs(criterion(2:end-1)) > threshold;
    increasing = abs(criterion(3:end)) < abs(criterion(2:end-1));
    decreasing = abs(criterion(1:end-2)) < abs(criterion(2:end-1));
    events = [false, thres & increasing & decreasing, false];
end

function scale = get_scale(template, data, N)
    scale_numerator = sum(template.*data) - sum(template).*sum(data)/N;
    scale_denominator = sum(template.^2) - sum(template).^2/N;
    scale = scale_numerator ./ scale_denominator;
end

function template = get_template(n, N, step, rise_time, decay_time)
    t = ((1-n):(N-n)) * step;
    template = t * 0;
    norm_factor = 1; % No idea what this is for and don't think it does anything anyway since 'scale' will account for whatever value this is.
    template(t > 0) = norm_factor * (1 - exp(-t(t > 0)/rise_time)) .* exp(-t(t > 0)/decay_time);
end
