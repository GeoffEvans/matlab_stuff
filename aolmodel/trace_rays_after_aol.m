function [zFocusModel] = trace_rays_after_aol(rayBundle, isPointingModeAndSingleBundle)

zFocusModel = FindModelFocus(rayBundle.zFocusExpected, isPointingModeAndSingleBundle);

normalToPlane = repmat([0 0 1]',1,rayBundle.numOfRays);
xyz = rayBundle.GetXyzLeavingAol();
k = rayBundle.k;

rayBundle.xyz{end-2} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(zFocusExpected,rayBundle));
rayBundle.xyz{end-1} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(zFocusModel,rayBundle));
rayBundle.xyz{end} = propagate_ray_to_plane(xyz,k,normalToPlane,PointAtZ(zFocusExpected*1.2,rayBundle));

    function point = PointAtZ(zVal,rayBundle)
        point = repmat([0 0 zVal]',1,rayBundle.numOfRays);
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