classdef aol_fft
    
    properties
        adjustment1 = 5e2; % ratio of k to max(kx)
        adjustment2 = 1e0;
        number_of_samples = 2.^8 - 1; % computational speed vs accuracy
        k = 2*pi/850e-9;
        z_list = linspace(60,90,99)*1e-6;
        do_plot = struct('input', 0, 'focal_plane', 0, 'angular_spec', 0);
        V = 613;
    end
    
    methods
        function [propagated_wave, x, y] = calculate_psf_through_aol(obj, waves, time)
            [sampled_wave_2d, space_width] = obj.get_sampled_wavefunction(waves, time);
            [propagated_wave, x, y] = obj.calculate_psf(sampled_wave_2d, space_width);
        end
        
        function [sampled_wave, space_width] = get_sampled_wavefunction(obj, waves, t)     
            num_aods = numel(waves.aods);
            space_width = pi * obj.number_of_samples / obj.k * obj.adjustment1;
            samples = linspace(-1/2, 1/2, obj.number_of_samples) * space_width;
            [x, y] = meshgrid(samples);

            distances = [0, ones(1,num_aods-1)*0.0]; %qq
            directions = linspace(0, 2 * pi, num_aods + 1);
            T = obj.V * t;
            r = arrayfun(@(theta) x.*cos(theta) + y.*sin(theta), directions, 'uniformoutput', 0);
            envelope = exp(-(x.^2 + y.^2)/2/(4e-3).^2);
            labels = {'x', 'y'};
            
            sampled_wave = envelope;
            for n = 1:num_aods
                phase = 2*pi * waves.aods{n} ./ (6e-3).^(1:5); % take aperture to be 12mm because of scaling by relay
                sampled_wave = obj.propagate_wave(sampled_wave, distances(n), space_width, 1)...
                    .* exp(1i * (phase(1) .* (r{n}-T) + phase(2) .* (r{n}-T).^2 + phase(3) .* (r{n}-T).^3 + phase(4) .* (r{n}-T).^4));
                if obj.do_plot.input
                    obj.plot_wavefunction_2d(sampled_wave, x, y, labels)
                end
            end

            focal_length = inv(num_aods * phase(2) / obj.k)
            sampled_wave = sampled_wave .* (sqrt(x.^2 + y.^2) < 6e-3) .* exp(1i * (...  
                (x.^2 + y.^2) * waves.focus * 2*pi ./ (6e-3).^2 ./ 2 ...  
                + tand(waves.ac_angle) * 0.5 .* (x.^3 + y.^3) .* focal_length.^-2  * obj.k/(2*pi)...  
                + (x.^2 + y.^2).^2 ./ (6e-3).^4 * waves.spherical * 2*pi)); 
            
            if obj.do_plot.input
                obj.plot_wavefunction_2d(sampled_wave, x, y, labels)
            end
        end
        
        function [propagated_wave, x, y] = calculate_psf(obj, sampled_wave_2d, space_width)
            [x, y, focal_plane_wave_2d, space_width] = obj.get_focal_plane_wavefunction(sampled_wave_2d, space_width);
            propagated_wave = obj.propagate_wave(focal_plane_wave_2d, obj.z_list, space_width, 4/3);
        end

        function [u, v, focal_plane_wave_2d, space_width] = get_focal_plane_wavefunction(obj, sampled_wave_2d, space_width)
            focal_plane_wave_2d = fftshift(fft2(ifftshift(sampled_wave_2d))); % transform to focal plane
            integer_list = floor(-obj.number_of_samples/2+0.5):floor(obj.number_of_samples/2-0.5);
            [u, v] = meshgrid(integer_list * 2 * pi / space_width / obj.k * 8e-3); % new positions in the focal plane   
            space_width = (max(u(:)) - min(u(:))) * obj.adjustment2;
            if obj.do_plot.focal_plane
                obj.plot_wavefunction_2d(focal_plane_wave_2d, u, v, {'x focus', 'y focus'})
            end
        end
        
        function propagated_wave = propagate_wave(obj, wave_2d, distance, space_width, ref_ind)
            count = [0:floor((obj.number_of_samples - 1)/2) floor(-(obj.number_of_samples - 1)/2):-1];
            [k_x, k_y] = meshgrid(count * 2 * pi / space_width);
            ft_wave_2d = fft2(ifftshift(wave_2d)); % recentre wave with fftshift

            %display(max(k_x(:))/obj.k); display((k_y(2)-k_y(1))/obj.k);
            k_z = sqrt(obj.k.^2 * ref_ind.^2 - k_x.^2 - k_y.^2);
            mask = ~(abs(imag(k_z)) > 0);

            propagated_wave = zeros([size(k_z), numel(distance)]);
            for n = 1:numel(distance)
                propagated_wave(:,:,n) = fftshift(ifft2((ft_wave_2d .* exp(1i * mask .* k_z .* distance(n)) .* mask))); % shift back 
            end
            propagated_wave = squeeze(propagated_wave);
            
            if obj.do_plot.angular_spec
                obj.plot_wavefunction_2d(fftshift(ft_wave_2d), fftshift(k_x), fftshift(k_y), {'kx', 'ky'})
            end
        end
        
        function plot_wavefunction_2d(~, wave_function_2d, x, y, labels)
            figure()
            subplot(1,2,1)
            hh = pcolor(x,y,abs(wave_function_2d));
            xlabel(labels{1})
            ylabel(labels{2})
            set(hh,'EdgeColor','none')

            subplot(1,2,2)
            hh = pcolor(x,y,angle(wave_function_2d));
            xlabel(labels{1})
            ylabel(labels{2})
            set(hh,'EdgeColor','none')
        end
        
        function xy_xz_plot_psf(obj, propagated_wave_2d, x, y, z_plane_frac)
            figure(); 
            subplot(1,2,1)
            h = pcolor(x, y, abs(propagated_wave_2d(:,:,round(z_plane_frac*size(obj.z_list,1)))).^1);
            set(h,'EdgeColor','none')
            axis equal

            subplot(1,2,2)
            [zz, xx] = meshgrid(obj.z_list, max(x,[],1));
            h = pcolor(xx, zz, abs(squeeze(propagated_wave_2d(:,round(size(x,2)/2),:))).^1);
            set(h,'EdgeColor','none')
        end
    end
end
