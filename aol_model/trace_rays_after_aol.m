function trace_rays_after_aol(rb, isPointingModeAndSingleBundle)

enlargementArray = [1,rb.numOfRaysPerPerturbation,rb.numOfPerturbations];
normalToPlane = repmat([0; 0; 1],enlargementArray);
xyz = rb.GetXyzLeavingAol();
k = rb.k;

[rb.zFocusModel, rb.xyFocusModel] = FindModelFocus(rb.zFocusPredicted, isPointingModeAndSingleBundle);

rb.xyz{end-2} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(rb.zFocusPredicted));
rb.xyz{end-1} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(rb.zFocusModel));
rb.xyz{end} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(rb.zFocusPredicted*1.2));

    function point = PointAtZ(zVal)
        point = repmat([0 0 zVal]',enlargementArray);
    end

    function [zFocus, xyFocus] = FindModelFocus(zFocusExpected, isPointingModeAndSingleBundle)
        zFocus = zFocusExpected;
        
        if isPointingModeAndSingleBundle == true
            zFocus = fminsearch(@MinFunc,zFocusExpected);
        end
        
        xyzFocus = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(zFocus));
        xyFocus = [mean(xyzFocus(1,:)), mean(xyzFocus(2,:))];
        
        function val = MinFunc(zVal)
            xyzTemp = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(zVal));
            sigmaX = std(xyzTemp(1,:),1);
            sigmaY = std(xyzTemp(2,:),1);
            val = prod(sigmaX) .* prod(sigmaY);
        end
    end
end