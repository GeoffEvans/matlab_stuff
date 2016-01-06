classdef TestCalcDriveForLvWrapper < matlab.unittest.TestCase
    
    properties
    end
    
    methods
        function tst = TestCalcDriveForLvWrapper()            
        end
    end
    
    methods(Test)
        function test_norm(tst)
            [scan_var, system_var, xyz_norm] = get_test_vars();
            [ a, b, c, d, ticks_per_ramp ] = calc_drive_for_lv_wrapper(scan_var, system_var, xyz_norm);
        end
    end
end