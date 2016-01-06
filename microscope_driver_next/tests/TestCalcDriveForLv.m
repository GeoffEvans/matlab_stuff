classdef TestCalcDriveForLv < matlab.unittest.TestCase
    
    properties
        aol
        scan
        obj_focal_len
        relay_mag
        num_rays
        num_rays_sqr
    end
    
    methods
        function tst = TestCalcDriveForLv()
            tst.aol = get_test_aol_params();
            
            s = get_test_scan_params();
            s.start_norm = [0.4; 0.2; 0.8];
            s.stop_norm = [0.4; 0.2; 0.8];
            tst.scan = s;
            
            tst.obj_focal_len = 8e-3;
            tst.relay_mag = 0.4;
            
            tst.num_rays = 3;
            tst.num_rays_sqr = tst.num_rays.^2;
        end
        
        function x = get_ray_points_3d(tst, drives, z_focus_obj, t)
            num_drives = size(drives, 2);
            
            k = zeros(2, tst.num_rays_sqr, num_drives);
            x = cell(7,1);
            x{1} = repmat(get_combinations(linspace(-1,1,tst.num_rays) * 1e-2), [1, 1, num_drives]);
            
            z_distance = [tst.aol.aod_reduced_z(2:end) - tst.aol.aod_reduced_z(1:(end-1)), 0];
            K_unit = tst.aol.aol_config.get_aod_dirs();
            aod_centre = tst.aol.aod_xy_offsets - repmat(tst.aol.aod1_centre_to_ref(1:2)', 1, 4);
            
            for m = 1:4
                T = t + tst.aol.transducer_centre_dist(m) ./ tst.aol.aod_ac_vel(m)...
                    - dot(x{m} - repmat(aod_centre(:,m), 1, tst.num_rays_sqr),...
                    repmat(K_unit(:,m), 1, tst.num_rays_sqr)) ./ tst.aol.aod_ac_vel(m);
                f = repmat(drives(1,:,m), 1, tst.num_rays_sqr) +...
                    repmat(drives(2,:,m), 1, tst.num_rays_sqr) .* T +...
                    repmat(drives(3,:,m), 1, tst.num_rays_sqr) .* T.^2 +...
                    repmat(drives(4,:,m), 1, tst.num_rays_sqr) .* T.^3;
                k = k - K_unit(:,m) * (tst.scan.op_wavelength / tst.aol.aod_ac_vel(m)) * f;
                x{m+1} = x{m} + k .* z_distance(m);
            end
            
            [x{5}, k] = tst.relay(x{5}, k);
            [x{6}, k] = tst.obj(x{5}, k);
            x{7} = x{6} + k .* z_focus_obj;
        end
        
        function [x_rel, k_rel] = relay(tst, x, k)
            mag = tst.relay_mag;
            x_rel = - x .* mag;
            k_rel = - k ./ mag;
        end
        
        function [x_obj, k_obj] = obj(tst, x, k)           
            r = sqrt(x(1,:).^2 + x(2,:).^2);
            theta = atan2(x(2,:), x(1,:));
            phi = atan2(k(2,:), k(1,:));
            angle_in = atan(sqrt(k(1,:).^2 + k(2,:).^2));

            r2 = tst.obj_focal_len .* angle_in;
            angle_out = -asin(r ./ tst.obj_focal_len);

            x_obj = [r2 .* cos(phi); r2 .* sin(phi)];
            k_obj = tst.convert_angles_to_wavevec(angle_out, theta);
            
            if ~isreal(k_obj)
                display('make relay mag smaller')
                for n = 1:size(k_obj, 2)
                    if ~isreal(k_obj(:,n))
                        k_obj(:,n) = zeros(2,1);
                        x_obj(:,n) = zeros(2,1);
                    end
                end
            end
        end
        
        function k = convert_angles_to_wavevec(tst, ang, theta)
            k(1,:) = cos(theta) .* sin(ang);
            k(2,:) = sin(theta) .* sin(ang);
        end
        
        function [x_img, y_img, z_img] = convert_norm_to_img(tst, x_norm, y_norm, z_norm)
            m = tst.relay_mag;
            f = tst.obj_focal_len;
            
            x_img = - 2 .* tst.scan.acceptance_angle .* f ./ m .* x_norm;
            y_img = - 2 .* tst.scan.acceptance_angle .* f ./ m .* y_norm;
            z_img = - 4 .* tst.scan.acceptance_angle ./ tst.aol.aod_aperture(end) .* (f ./ m).^2 .* z_norm;
        end
    end
    
    methods(Test)
        function test_norm(tst)
            norms = tst.scan.start_norm;
            [x_img, y_img, z_img] = convert_norm_to_img(tst, norms(1,:), norms(2,:), norms(3,:));
            
            drives = calc_drive_for_lv(tst.aol, tst.scan);
            expected = repmat([x_img;y_img], 1, tst.num_rays_sqr);
            
            x1 = get_ray_points_3d(tst, drives, z_img, 0);
            tst.assertEqual(x1{end}, expected, 'RelTol', 10^-2);
            
            x2 = get_ray_points_3d(tst, drives, z_img, -2e-6);
            tst.assertEqual(x2{end}, expected, 'RelTol', 10^-2);
        end
    end
end