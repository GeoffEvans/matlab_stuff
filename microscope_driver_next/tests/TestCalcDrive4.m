classdef TestCalcDrive4 < matlab.unittest.TestCase
    
    properties
        aol
        scan
        num_rays = 5;
        num_rays_sqr
    end
    
    methods(Test)
        function test_simple_scan(tst)
            x_pos = 0.02;
            start = [x_pos;x_pos;1];
            stop = start;
            drive_coeffs_aod_centres = calc_drive_4(tst.aol, start, stop, tst.scan);
            
            x1 = tst.get_ray_points_3d(drive_coeffs_aod_centres, 1, 0);
            tst.assertEqual(x1{end}, repmat(x_pos, size(x1{end})), 'RelTol', 1e-2);
            
            x2 = tst.get_ray_points_3d(drive_coeffs_aod_centres, 1, 1e-6);
            tst.assertEqual(x2{end}, repmat(x_pos, size(x2{end})), 'RelTol', 1e-2);
        end
    end
    
    methods
        function t = TestCalcDrive4()
            t.aol = get_test_aol_params();
            t.scan = get_test_scan_params();
            t.num_rays_sqr = t.num_rays.^2;
        end
        
        function [x] = get_ray_points_3d(tst, drives, z_focus, t)
            num_drives = size(drives, 2);
            
            k = zeros(2, tst.num_rays_sqr, num_drives);
            x = cell(5,1);
            x{1} = repmat(get_combinations(linspace(-1,1,tst.num_rays) * 1e-2), [1, 1, num_drives]);
            
            z_distance = [tst.aol.aod_reduced_z(2:end) - tst.aol.aod_reduced_z(1:(end-1)), z_focus];
            K_unit = tst.aol.aol_config.get_aod_dirs();
            aod_centre = tst.aol.aod_xy_offsets - repmat(tst.aol.aod1_centre_to_ref(1:2)', 1, 4);
            
            for m = 1:4
                T = t - dot(x{m} - repmat(aod_centre(:,m), 1, tst.num_rays_sqr),...
                    repmat(K_unit(:,m), 1, tst.num_rays_sqr)) ./ tst.aol.aod_ac_vel(m);
                f = repmat(drives(1,:,m), 1, tst.num_rays_sqr) +...
                    repmat(drives(2,:,m), 1, tst.num_rays_sqr) .* T +...
                    repmat(drives(3,:,m), 1, tst.num_rays_sqr) .* T.^2 +...
                    repmat(drives(4,:,m), 1, tst.num_rays_sqr) .* T.^3;
                k = k - K_unit(:,m) * (tst.scan.op_wavelength / tst.aol.aod_ac_vel(m)) * f;
                x{m+1} = x{m} + k .* z_distance(m);
            end
        end
    end
end