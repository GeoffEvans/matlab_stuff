function JC2()

wavefront_with_rays()
end

function wavefront_with_rays()
	x = linspace(-1,1,1000);
    x_sparse = linspace(-0.9, 0.9, 10);
    z1 = @(x) - x/5 + x.^2/10 + x.^3/8 - x.^4/8;
    n1 = @(x) - 1/5 + x/5 + x.^2/8*3 - x.^3/2;
    scale = @(x) 1./sqrt(1+n1(x).^2);
    
    hold on
    plot(x,z1(x), 'r')
    plot([x_sparse; x_sparse - 1 * scale(x_sparse) .* n1(x_sparse)], [z1(x_sparse); z1(x_sparse) + 1 * scale(x_sparse)], 'b')
    plot(x - 1 * scale(x) .* n1(x), z1(x) + 1 * scale(x), 'r');
    plot([x_sparse; x_sparse - 2 * scale(x_sparse) .* n1(x_sparse)], [z1(x_sparse); z1(x_sparse) + 2 * scale(x_sparse)], 'b')
    plot(x - 2 * scale(x) .* n1(x), z1(x) + 2 * scale(x), 'r');
    hold off
    axis([-2,1.5,-0.5,2.5])
    axis equal
end

function quad_vs_circ(plot_quad_rays, plot_circ_rays)
    x = linspace(-1,1,1000);
    x_sparse = linspace(-0.6, 0.6, 8);
    z1 = @(x) x.^2/2;
    n1 = @(x) x;
    z2 = @(x) 1 - sqrt(1 - x.^2);
    n2 = @(x) x ./ sqrt(1 - x.^2);
    
    hold on
    plot(x,z1(x), 'k')
    plot(x,z2(x), 'r')
    
    if plot_quad_rays
        plot([x_sparse; x_sparse - 1.5 * n1(x_sparse)], [z1(x_sparse); z1(x_sparse) + 1.5], 'b')
    end
    if plot_circ_rays
        plot([x_sparse; x_sparse - 1.5 * n2(x_sparse)], [z2(x_sparse); z2(x_sparse) + 1.5], 'b')
    end
    
    axis equal
    hold off
end

