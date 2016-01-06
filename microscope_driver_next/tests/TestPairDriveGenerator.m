classdef TestPairDriveGenerator < matlab.unittest.TestCase
    
    properties
        p
        num_rays = 5;
    end
    
    methods(Test)
        function test_quad_stationary(t)
            D = 0;
            z = [1, 2];
            x = [1e-2, -1e-2];
            
            t.do_test(x, x, z, z, D, 3)
        end
        function test_quad_moving(t)
            D = 0;
            z = [1, 2];
            x_start = [-1e-2, -1e-2];
            x_end = [1e-2, -5e-3];
            
            t.do_test(x_start, x_end, z, z, D, 2);
        end
        function test_cubic_z(t)
            D = 0;
            x = [1e-2, -1e-2];
            z_start = [1, 2];
            z_end = [1.2, 1.6];
            
            t.do_test(x, x, z_start, z_end, D, 3)      
        end
        function test_cubic_xz(t)
            D = 0;
            x_start = [1e-2, -1e-2];
            x_end = [-1e-2, -5e-3];
            z_start = [1, 2];
            z_end = [1.2, 1.6];
            
            t.do_test(x_start, x_end, z_start, z_end, D, 2)      
        end
        function test_quartic_simple(t)
            D = 1e2;
            z = [1, 2];
            x = [1e-2, -1e-2];
            
            t.do_test(x, x, z, z, D, 1.5)       
        end
        function test_quartic(t)
            D = 1e1;
            x_start = [1e-2, -1e-2];
            x_end = [-1e-2, -5e-3];
            z_start = [1, 2];
            z_end = [1.2, 1.6];
            
            t.do_test(x_start, x_end, z_start, z_end, D, 2)
        end
    end
    
    methods 
        function do_test(t, x_start, x_end, z_start, z_end, D, n)
            t.create_pair(x_start, x_end, z_start, z_end, D);
            [drive1, drive2] = t.p.get_coefficients_at_ref();
            [~, ~, ray_x_start] = t.get_ray_points(drive1, drive2, z_start, -t.p.time_half);
            [~, ~, ray_x_end] = t.get_ray_points(drive1, drive2, z_end, t.p.time_half);
            
            %simulate_rays(t, drive1, drive2, z_end); % for graphically checking test 
            
            t.assertEqual(repmat(x_start,5,1), ray_x_start, 'RelTol', 0.1^n);
            t.assertEqual(repmat(x_end,5,1), ray_x_end, 'RelTol', 0.1^n);
        end
        
        function [x1, x2, x3] = get_ray_points(test, drive1, drive2, z, t)
            num_drives = size(drive1, 2);
            
            k = zeros(test.num_rays, num_drives);
            x1 = repmat(linspace(-1,1,test.num_rays)', 1, num_drives) * 1e-2;

            T = t - x1 / test.p.V;
            f1 = repmat(drive1(1,:), test.num_rays, 1) +...
                repmat(drive1(2,:), test.num_rays, 1) .* T +...
                repmat(drive1(3,:), test.num_rays, 1) .* T.^2 +...
                repmat(drive1(4,:), test.num_rays, 1) .* T.^3;
            k = k - (test.p.lambda / test.p.V) * f1;
            x2 = x1 + k .* test.p.L;

            T = t + x2 / test.p.V;
            f2 = repmat(drive2(1,:), test.num_rays, 1) +...
                repmat(drive2(2,:), test.num_rays, 1) .* T +...
                repmat(drive2(3,:), test.num_rays, 1) .* T.^2 +...
                repmat(drive2(4,:), test.num_rays, 1) .* T.^3;
            k = k + (test.p.lambda / test.p.V) .* f2;
            x3 = x2 + k .* repmat(z, test.num_rays, 1);
        end
        
        function simulate_rays(tst, drive1, drive2, z)
            figure();
            for n = -1:1
                t = n * tst.p.time_half;
                
                [x1, x2, x3] = tst.get_ray_points(drive1, drive2, z, t);
                
                hold on;
                for col = 1:size(x1, 2)
                    z_pos = cumsum([0, 0.1, tst.p.L, z(col)]);
                    for u = 1:size(x1, 1)
                        plot(z_pos,[x1(u,col),x1(u,col),x2(u,col),x3(u,col)], 'color', (n + 2)/5 * [1, -1, 0] + [0, 1, 0]);
                    end
                end
                hold off;
            end
        end
        
        function create_pair(t, x_start, x_stop, z_start, z_stop, D)
            t.p = PairDriveGenerator(800e-9, 613, 5e-2, 40e6,...
                x_start, x_stop, z_start, z_stop, 1e-5, D);
        end
    end
end