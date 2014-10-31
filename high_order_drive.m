function high_order_drive()

w = [0, 0, 0, 0, 0];
w(1) = 0; % angle of deflection
w(2) = 1; % 1/focal length
w(3) = 0; % coma -2
w(4) = 1; % spherical 0
w(5) = 0; % spherical 3

c = [0, 1e11, 1e15, 0, 0];
w = wavefront_to_chirps(c, false, 1);
w(4) = -140;
w = propagate_wavefront(w, 0.05);
w(3) = 0;

plot_wavefront(w)

end

function plot_wavefront(w)
    aod_width = 20e-3;
    x_array = linspace(-aod_width/2, aod_width/2, 20);
    z_array = wavefront_func(x_array, w);

    subplot(1,3,1)
    plot(z_array, x_array)
    axis tight
    xlabel('z')
    ylabel('x')
    
    subplot(1,3,2:3)
    add_rays(x_array, w)
end

function add_rays(x, w)
    cla
    hold on
    focal_length = 1/w(2);    

    ray_origin_z = wavefront_func(x, w);
    ray_origin_x = x; 

    shift = 1e-4;
    ray_start_z = ray_origin_z + focal_length - shift;
    ray_start_x = ray_origin_x + (focal_length - shift) .* wavefront_normal(x, w);
    ray_end_z = ray_origin_z + focal_length + shift;
    ray_end_x = ray_origin_x + (focal_length + shift) .* wavefront_normal(x, w);

    plot([ray_start_z; ray_end_z], [ray_start_x; ray_end_x])
    axis tight
    hold off
end

function z = wavefront_func(x, w)
    z = w(1)*x + w(2)/2*x.^2 + w(3)/6*x.^3 + w(4)/24*x.^4 + w(5)/120*x.^5;
end

function normal_grad = wavefront_normal(x, w)
    dzdx = w(1) + w(2)*x + w(3)/2*x.^2 + w(4)/6*x.^3 + w(5)/24*x.^4;
    normal_grad = -1 * dzdx;
end

function w_new = propagate_wavefront(w, d)
    B = 1/(1-d*w(2));
    w_new(1) = w(1);
    w_new(2) = B * w(2);
    w_new(3) = B^3 * w(3);
    w_new(4) = B^4 * ( w(4) + 3*d*(B*w(3)^2 - w(2)^4) );
    w_new(5) = B^5 * ( w(5) + 5*B*d*w(3) * (2*w(4) + 3*B*d*w(3)^2 - 6*w(2)^3) );
end

function out = wavefront_to_chirps(in, w_to_c, dir)
    V = 613;
    lambda = 800e-9;
    
    if w_to_c
        w = in;
        out = [-dir*V*w(1), V^2*w(2), -dir*V^3*w(3)/2, V^4*w(4)/6, -dir*V^5*w(5)/24] / lambda;
    else
        c = in;
        out = [-dir/V*c(1), 1/V^2*c(2), -dir/V^3*c(3)*2, 1/V^4*c(4)*6, -dir/V^5*c(5)*24] * lambda;
    end
end