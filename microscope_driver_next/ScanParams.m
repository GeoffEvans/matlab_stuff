classdef ScanParams
    properties
        start_norm
        stop_norm
        ramp_time
        
        aod_opt_freq
        acceptance_angle
        op_wavelength
        
        zoom_factor
        voxel_density_1D
        waves_of_4th
    end
    methods
        function s = ScanParams(scan_var, xyz_norm)
            if nargin == 0
                return
            end
            
            s.op_wavelength = scan_var.opticalWaveLength * 1e-9; % convert from nm to m
            s.acceptance_angle = scan_var.acceptanceAngle;            
            s.aod_opt_freq = repmat(scan_var.centreFrequency * 1e6, 1, 4); % convert to MHz to Hz   
            s.zoom_factor = scan_var.zoomFactor;
            s.voxel_density_1D = scan_var.numberOfVoxels;
            display('waves of 4th order aberration hardcoded')
            s.waves_of_4th = [0, 0]; % TODO
            
            aol_modes = [ImagingMode.Raster, ImagingMode.Structural, ImagingMode.Functional, ImagingMode.Miniscan];
            imaging_mode = aol_modes(scan_var.aolMode);
            [s.start_norm, s.stop_norm] = s.handle_image_start_stops(xyz_norm, imaging_mode);
            
            pixel_time = scan_var.dwellTime;
            s.ramp_time = s.calc_ramp_time(pixel_time);
        end
        
        function D = get_D(scan, aol, idx)
            D = scan.op_wavelength .* scan.waves_of_4th(idx) .* (aol.aod_aperture(end)/2).^-4;
        end
        
        function ramp_time = calc_ramp_time(s, pixel_time)
            ramp_time = pixel_time .* s.voxel_density_1D .* mag(s.start_norm - s.stop_norm) ./ 2; % time per pixel * num pixels
        end
        
        function [start, stop] = handle_image_start_stops(s, xyz_norm, imaging_mode)
            if imaging_mode == ImagingMode.Raster
                [start, stop] = generate_rows(s, xyz_norm);
            elseif imaging_mode == ImagingMode.Structural
                [start, stop] = generate_grid(s, xyz_norm);
            else
                start = xyz_norm.imageStartNormalised';
                stop = xyz_norm.imageStopNormalised';
            end
        end

        function [start, stop] = generate_grid(s, xyz_norm)
            start_norm_in = xyz_norm.imageStartNormalised;
            line_pts = linspace(-1, 1, s.voxel_density_1D) ./ s.zoom_factor;
            start_x = start_norm_in(1) + line_pts;
            start_y = start_norm_in(2) + line_pts;
            start_z = start_norm_in(3) * ones(size(start_x));
            start = [start_x; start_y; start_z];
            stop = start;
        end

        function [start, stop] = generate_rows(s, xyz_norm)
            start_norm_in = xyz_norm.imageStartNormalised;
            line_y = linspace(-1, 1, s.voxel_density_1D) ./ s.zoom_factor;
            start_x = start_norm_in(1) - ones(size(line_y)) ./ s.zoom_factor;
            stop_x = start_norm_in(1) + ones(size(line_y)) ./ s.zoom_factor;
            start_y = start_norm_in(2) + line_y;
            start_z = start_norm_in(3) * ones(size(start_x));
            start = [start_x; start_y; start_z];
            stop = [stop_x; start_y; start_z];
        end
    end
end
    

