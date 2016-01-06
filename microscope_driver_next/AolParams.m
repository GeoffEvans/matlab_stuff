classdef AolParams
    properties
        num_aods
        aol_config
        system_clock_freq
        data_time_interval
        
        aod_ac_vel
        aod_aperture
        transducer_centre_dist

        aod_mode
        aod_xy_offsets
        aod1_centre_to_ref
        aod_reduced_z
    end
    
    methods
        function d = AolParams(system_var)
            if nargin == 0
                return
            end
            
            d.aod_ac_vel = repmat(system_var.acousticVelocity, 1, 4);
            d.aod_mode = system_var.diffractionMode;
            d.aod_aperture = system_var.crystalLength;  
            d.system_clock_freq = system_var.controlClockFreq;
            d.data_time_interval = system_var.dataTimeInterval;

            disp('hardcoded aol config - TODO')
            d.aol_config = AolConfig.Orig_4;
            d.num_aods = d.aol_config.get_num_aods();
            disp('hardcoded aod xy offsets - TODO')
            d.aod_xy_offsets = [0, 0.0026, 0.0052, 0.0052; 0, 0, 0.0026, 0.0052];
            disp('hardcoded centre reference displacement - TODO')
            d.aod1_centre_to_ref = [-0.005, -0.005, 16e-2];
            disp('hardcoded transducer centre distances - TODO')
            d.transducer_centre_dist = [8, 8, 8, 8] * 1e-3;           

            disp('hardcoded ord ref index - TODO')
            ord_ref_ind = 2.26;
            aod_thickness = [system_var.thicknessOfAod1, system_var.thicknessOfAod2, system_var.thicknessOfAod3, system_var.thicknessOfAod4];
            aod_spacing = [system_var.seperationOfAod1Aod2, system_var.seperationOfAod2Aod3, system_var.seperationOfAod3Aod4];
            d.aod_reduced_z = d.calc_reduced_z(aod_spacing, ord_ref_ind, aod_thickness);            
        end
                    
        function aod_reduced_z = calc_reduced_z(d, aod_spacing, ord_ref_ind, aod_thickness)
            aod_z = [0, cumsum(aod_spacing)] - d.aod1_centre_to_ref(3);
            thickness_reduction = (1 - ord_ref_ind.^-1) .* aod_thickness;
            aod_reduced_z = aod_z + fliplr(cumsum(fliplr(thickness_reduction)));
        end
        
        function [a, b, c, d, ticks_per_ramp] = get_params_for_lv(params, drive_coeffs_transducers, ramp_time)    
            swap = params.aol_config == AolConfig.Orig_4;
            
            aod_fill_time = params.data_time_interval .* max(ceil(params.aod_aperture ./ params.aod_ac_vel ./ params.data_time_interval));
            scan_time = aod_fill_time + ramp_time;
            
            ticks_per_ramp = round(scan_time .* params.system_clock_freq);
            a = squeeze_swap_scale_round(drive_coeffs_transducers(1,:,:), 2^32 / params.system_clock_freq);
            b = squeeze_swap_scale_round(drive_coeffs_transducers(2,:,:), 2^32 / params.system_clock_freq^2 / 2^3); % scale B down here and expand later
            c = squeeze_swap_scale_round(drive_coeffs_transducers(3,:,:), 2^32 / params.system_clock_freq^3 * 2^13); 
            d = squeeze_swap_scale_round(drive_coeffs_transducers(4,:,:), 2^32 / params.system_clock_freq^4 * 2^25);
            display('TODO: should these coefficients be scaled to be phases for the fpga?')
            
            function y = squeeze_swap_scale_round(x, scale)
                y = round(squeeze(x) * scale);
                if swap
                    y = y(:,[1 3 2 4]);
                end
            end
        end
    end
end

