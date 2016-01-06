classdef TestScanParams < matlab.unittest.TestCase
    
    properties
        scan
        start_norm
        stop_norm
        xyz_norm
    end
    
    methods(Test)
        function test_structural(tst)
            [start, stop] = tst.scan.handle_image_start_stops(tst.xyz_norm, ImagingMode.Structural);
            
            tst.assertEqual(start, repmat(tst.start_norm', 1, 3) + [-1, 0, 1; -1, 0, 1; 0, 0, 0]);
            tst.assertEqual(stop, repmat(tst.start_norm', 1, 3) + [-1, 0, 1; -1, 0, 1; 0, 0, 0]);
        end
        function test_raster(tst)
            [start, stop] = tst.scan.handle_image_start_stops(tst.xyz_norm, ImagingMode.Raster);
            
            tst.assertEqual(start, repmat(tst.start_norm', 1, 3) + [-1, -1, -1; -1, 0, 1; 0, 0, 0]);
            tst.assertEqual(stop, repmat(tst.start_norm', 1, 3) + [1, 1, 1; -1, 0, 1; 0, 0, 0]);
        end
        function test_miniscan(tst)
            [start, stop] = tst.scan.handle_image_start_stops(tst.xyz_norm, ImagingMode.Functional);
            
            tst.assertEqual(start, tst.start_norm');
            tst.assertEqual(stop, tst.start_norm');
        end
        function test_functional(tst)
            [start, stop] = tst.scan.handle_image_start_stops(tst.xyz_norm, ImagingMode.Functional);
            
            tst.assertEqual(start, tst.start_norm');
            tst.assertEqual(stop, tst.start_norm');
            tst.assertNotEqual(stop, tst.start_norm);
        end
    end
    
    methods 
        function tst = TestScanParams()
            s = get_test_scan_params();
            s.voxel_density_1D = 3;
            tst.scan = s;
            
            tst.start_norm = [0, 0.1, 0.5];
            tst.stop_norm = [0, 0.1, 0.5];
            tst.xyz_norm = struct('imageStartNormalised', tst.start_norm, 'imageStopNormalised', tst.stop_norm);
        end
    end
end