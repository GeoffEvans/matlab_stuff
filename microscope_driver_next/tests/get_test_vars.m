function [scan_var, system_var, xyz_norm] = get_test_vars()

display('make sure the params are in cell arrays')
scan_raw = [1.0000 800 40 4.3500 1.0000 200.0000 0.0800];
system_raw = [0 0 0 0 240e6 200e6 50e-9 -1 15 613 0.05 0.05 0.05 0.008 0.008 0.08 0.08 1];
xyzImageStartNormalised = [0 0 0];
xyzImageStopNormalised = [0 0 0];

[scan_var, xyz_norm, system_var] = variableConstruct(scan_raw, xyzImageStartNormalised, xyzImageStopNormalised, system_raw);
end

