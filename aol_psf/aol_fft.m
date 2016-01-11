classdef aol_fft
    
    properties
        adjustment = 1e2; % adjustments as necessary for the FFT - sets ratio of k to max(kx)
        number_of_samples = 2.^8 - 1; % computational speed vs accuracy
        k = 2*pi/800e-9; % wavevector
        z_list = linspace(-100,100,199)*1e-6; % distances from the nominal focal plane of the objective to evaluate the field at 
        do_plot = struct('aods', 0, 'input', 0, 'focal_plane', 0, 'angular_spec', 0); % plotting of the intermediate stages useful for debugging
        V = 613; % speed of sound in TeO2
        half_width = 7.5e-3; % half the aperture width
        scaling = 0.8; % scaling by the relay between AOL and objective
        mag = 20;
        tube_focal_len = 160e-3;
        beam_sigma = 3e-3; % 68% of field within +- beam_sigma, 95% of field within +- 2 beam_sigma, equiv. beam intensity falls off to 1/e after beam_sigma
        pwr = 1;
        cm = 'jet';
        spacing = 0.04554;
        ramp_width = 1e9;
        number_of_ramps = 1;
    end
    
    methods
        function res = get_psf(obj, time, ws, wf, aod_drives_in_waves, plot)
            waves = struct();
            waves.aods = aod_drives_in_waves;
            waves.focus = wf;
            waves.spherical = ws;
            waves.ac_angle = 0; % walk off angle in degrees

            [propagated_wave, x, y] = obj.calculate_psf_through_aol(waves, time);
            if plot
                figure()
                subplot(1,2,1); obj.plot_psf_xy(propagated_wave, x, y);
                subplot(1,2,2); obj.plot_psf_xz(propagated_wave, x);
            end
            res = obj.get_psf_dimensions(propagated_wave, x, y);
        end

        function [propagated_wave_max, x, y] = calculate_psf_through_aol(obj, waves, time)
            % take a number of waves of phase shift for each AOD and a time to calculate the PSF
            propagated_wave_max = 0;
            if obj.k == 2*pi / 800e-9
                wavelen_list = [800, 800, 800-2.5, 800+2.5] * 1e-9;
            elseif obj.k == 2*pi / 920e-9
                wavelen_list = [920, 920, 920-3.5, 920+3.5] * 1e-9;
            else
                wavelen_list = 2*pi / obj.k;
            end
            for wave_len = wavelen_list
                obj.k = 2*pi/wave_len;
                for t = time; 
                    [sampled_wave_2d, space_width] = obj.get_sampled_wavefunction(waves, t);
                    [propagated_wave, x, y] = obj.calculate_psf(sampled_wave_2d, space_width);
                    propagated_wave = propagated_wave ./ max(propagated_wave(:));
                    propagated_wave_max = max(propagated_wave_max, propagated_wave);
                end
            end
        end
        
        function ramp = get_ramp(obj, T, phase, r, n)
            freq = exp(1i * (phase(1) .* (r{n}-T) + phase(2) .* (r{n}-T).^2 + phase(3) .* (r{n}-T).^3 + phase(4) .* (r{n}-T).^4 + phase(5) .* (r{n}-T).^5));
            amp = exp(-(r{n}-T).^2/obj.ramp_width.^2);
            ramp = amp .* freq;
        end
            
        function [sampled_wave, space_width] = get_sampled_wavefunction(obj, waves, t)     
            % take a number of waves of phase shift for each AOD and a time to calculate the wave out of the last AOD
            num_aods = numel(waves.aods);
            space_width = pi * obj.number_of_samples / obj.k * obj.adjustment;
            samples = linspace(-1/2, 1/2, obj.number_of_samples) * space_width;
            [x, y] = meshgrid(samples);

            distances = [0, ones(1,num_aods-1)*obj.spacing];
            directions = linspace(0, 2 * pi, num_aods + 1);
            T = obj.V * t;
            r = arrayfun(@(theta) x.*cos(theta) + y.*sin(theta), directions, 'uniformoutput', 0);
            gaussian = exp(-(x.^2 + y.^2)/2/(obj.beam_sigma).^2);
            aperture_in = (sqrt(x.^2 + y.^2) < obj.half_width * 1.2);
            aperture_out = (sqrt(x.^2 + y.^2) < obj.half_width);
            labels = {'n', 'x', 'y'};
            
            sampled_wave = gaussian .* aperture_in;
            for n = 1:num_aods
                phase = 2*pi * waves.aods{n} ./ obj.half_width .^ (1:5);
                phase_shift = 0;
                for m = (1:obj.number_of_ramps) - 1/2 - obj.number_of_ramps/2 
                    phase_shift = phase_shift + obj.get_ramp(T - m*obj.ramp_width, phase, r, n);
                end
                sampled_wave = obj.propagate_wave(sampled_wave, distances(n), space_width, 1) .* phase_shift;
                
                if obj.do_plot.aods
                    labels{1} = n;
                    obj.plot_wavefunction_2d(sampled_wave, x, y, labels)
                end
            end
            sampled_wave = sampled_wave .* aperture_out .* exp(1i * (...  
                (x.^2 + y.^2) * waves.focus * 2*pi ./ obj.half_width .^ 2 ...  
                + (x.^2 + y.^2).^2 ./ obj.half_width .^4 * waves.spherical * 2*pi)); 
            
            if obj.do_plot.input
                figure; hold; plot(unwrap(ifftshift(angle(sampled_wave(ceil(size(sampled_wave,2)/2),:)))), 'k'); plot(unwrap(ifftshift(angle(sampled_wave(:,ceil(size(sampled_wave,2)/2))))), 'r--');
                labels{1} = 'input to relay and obj';
                obj.plot_wavefunction_2d(sampled_wave, x, y, labels)
            end
        end
        
        function [propagated_wave, x, y] = calculate_psf(obj, sampled_wave_2d, space_width)
            % taking the field into the objective, calculate the field over a volume to find the PSF
            [x, y, focal_plane_wave_2d, space_width] = obj.get_focal_plane_wavefunction(sampled_wave_2d, space_width);
            propagated_wave = obj.propagate_wave(focal_plane_wave_2d, obj.z_list, space_width, 4/3);
        end

        function [u, v, focal_plane_wave_2d, space_width] = get_focal_plane_wavefunction(obj, sampled_wave_2d, space_width)
            % FFT the wave into the (infinity corrected) objective to find the field in the nominal focal plane
            focal_plane_wave_2d = fftshift(fft2(ifftshift(sampled_wave_2d))); % transform to focal plane
            integer_list = floor(-obj.number_of_samples/2+0.5):floor(obj.number_of_samples/2-0.5);
            k_x_discrete = integer_list * 2 * pi / (space_width * obj.scaling);
            plane_wave_angle_air = k_x_discrete / obj.k;
            [u, v] = meshgrid(plane_wave_angle_air * obj.tube_focal_len / obj.mag); % new positions in the focal plane   
            space_width = (max(u(:)) - min(u(:)));
            if obj.do_plot.focal_plane
                obj.plot_wavefunction_2d(focal_plane_wave_2d, u, v, {'focal plane', 'x focus', 'y focus'})
            end
        end
        
        function propagated_wave = propagate_wave(obj, wave_2d, distance, space_width, ref_ind)
            % use FFT to propagate the field backwards and forwards
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
                obj.plot_wavefunction_2d(fftshift(ft_wave_2d), fftshift(k_x), fftshift(k_y), {'ang spec', 'kx', 'ky'})
            end
        end
        
        function plot_wavefunction_2d(~, wave_function_2d, x, y, labels)
            figure()
            subplot(1,2,1)
            hh = pcolor(x,y,abs(wave_function_2d));
            title(labels{1})
            xlabel(labels{2})
            ylabel(labels{3})
            set(hh,'EdgeColor','none')

            subplot(1,2,2)
            hh = pcolor(x,y,angle(wave_function_2d));
            title(labels{1})
            xlabel(labels{2})
            ylabel(labels{3})
            set(hh,'EdgeColor','none')
        end
        
        function plot_psf(obj, a, b, c)
            h = pcolor(a, b, c);
            set(h,'EdgeColor','none')
            colormap(obj.cm)
            %caxis([0,2000^obj.pwr])
            axis equal
            %axis tight
            %axis off
            %set(gcf, 'Position', [0,0,800,800]);
            %set(gca,'position',[0 0 1 1],'units','normalized')
        end
        function plot_psf_xy(obj, propagated_wave_2d, x, y)
            [~,z_plane_index] = max(squeeze(max(max(propagated_wave_2d))));
            obj.plot_psf(x, y, abs(propagated_wave_2d(:,:,z_plane_index)).^obj.pwr);
        end          
        function plot_psf_xz(obj, propagated_wave_2d, x)
            [zz, xx] = meshgrid(obj.z_list, max(x,[],1));
            obj.plot_psf(xx, zz, abs(squeeze(max(propagated_wave_2d.^obj.pwr, [], 1))));
        end            
        function plot_psf_yz(obj, propagated_wave_2d, y)
            [zz, yy] = meshgrid(obj.z_list, max(y,[],2));
            obj.plot_psf(yy, zz, abs(squeeze(max(propagated_wave_2d.^obj.pwr, [], 2))));
        end 
        
        function res = get_psf_dimensions(obj, field_3d, x, y)
            z = obj.z_list;
            % use FWHM measurements to quantify the PSF dimensions
            r = sqrt(x.^2 + y.^2);
            max_intensity_sqr = max(abs(field_3d(:)).^4);

            half_or_more_r = max(abs(field_3d), [], 3).^4 >= max_intensity_sqr/2;
            r_res = max(r(half_or_more_r)) - min(r(half_or_more_r));
            half_or_more_z = squeeze(max(max(abs(field_3d)))).^4 >= max_intensity_sqr/2;
            z_res = max(z(half_or_more_z)) - min(z(half_or_more_z));
            
            max_val = max(abs(field_3d(:)));
            r_pos = median(x(max(abs(field_3d), [], 3) == max_val));
            z_pos = median(z(max(max(abs(field_3d), [], 1), [], 2) == max_val));

            max_intensity_sqr = max(max(abs(field_3d(:,:,ceil(numel(z)/2))).^4));
            total_fluores = sum(abs(field_3d(:)).^4);
            res = [[r_pos, z_pos] * 1e6, [r_res, z_res] * 1e6, log10(max_intensity_sqr), log10(total_fluores)];
        end
    end
end
