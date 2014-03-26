function trace_rays_after_aol(rb, isPointingModeAndSingleBundle)

enlargementArray = [1,rb.numOfRaysPerPerturbation,rb.numOfPerturbations];
normalToPlane = repmat([0; 0; 1],enlargementArray);
xyz = rb.GetXyzLeavingAol();
k = rb.k;

rb.zFocusModel = FindModelFocus(rb.zFocusPredicted, isPointingModeAndSingleBundle);
fractionalFocusErrorZ = rb.zFocusPredicted/rb.zFocusModel - 1; % output a measure of focus accuracy

rb.xyz{end-2} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(rb.zFocusPredicted));
rb.xyz{end-1} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(rb.zFocusModel));
rb.xyz{end} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(rb.zFocusPredicted*1.2));

    function point = PointAtZ(zVal)
        point = repmat([0 0 zVal]',enlargementArray);
    end

    function focus = FindModelFocus(zFocusExpected, isPointingModeAndSingleBundle)
        focus = zFocusExpected;
        
        if isPointingModeAndSingleBundle == true
            focus = fminsearch(@MinFunc,zFocusExpected);
        end
        
        function val = MinFunc(zVal)
            xyzTemp = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(zVal));
            sigmaX = std(xyzTemp(1,:),1);
            sigmaY = std(xyzTemp(2,:),1);
            val = prod(sigmaX) .* prod(sigmaY);
        end
    end
end