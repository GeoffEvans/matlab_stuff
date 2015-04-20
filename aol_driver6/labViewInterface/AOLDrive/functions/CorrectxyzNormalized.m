 function [ xyzNormalised,xyzMid,xyzLength ] = ...
     CorrectxyzNormalized( xyzNormalised,deltaxyz,xyzLength)
% see whether for loop can be avoided.
numPoints = size(xyzNormalised.imageStartNormalised);
for i = 1:numPoints;
    if (xyzNormalised.imageStartNormalised(i,1) > ...
            xyzNormalised.imageStopNormalised(i,1))
        
        xyzNormalised.imageStartNormalised(i) = ...
            xyzNormalised.imageStartNormalised(i,1) + deltaxyz(i,1);
        
        xyzNormalised.imageStopNormalised(i,1) = ...
            xyzNormalised.imageStopNormalised(i,1) - deltaxyz(i,1);
    else
        xyzNormalised.imageStartNormalised(i,1) = ...
            xyzNormalised.imageStartNormalised(i,1) - deltaxyz(i,1);
        
        xyzNormalised.imageStopNormalised(i,1) = ...
            xyzNormalised.imageStopNormalised(i,1) + deltaxyz(i,1);
    end

    if (xyzNormalised.imageStartNormalised(i,2) > ...
            xyzNormalised.imageStopNormalised(i,2)) 
        
        xyzNormalised.imageStartNormalised(i,2) = ...
            xyzNormalised.imageStartNormalised(i,2) + deltaxyz(i,2);
        
        xyzNormalised.imageStopNormalised(i,2) = ...
            xyzNormalised.imageStopNormalised(i,2) - deltaxyz(i,2);
    else
        xyzNormalised.imageStartNormalised(i,2) = ...
            xyzNormalised.imageStartNormalised(i,2) - deltaxyz(i,2);
        
        xyzNormalised.imageStopNormalised(i,2) = ...
            xyzNormalised.imageStopNormalised(i,2) + deltaxyz(i,2);
    end

    if (xyzNormalised.imageStartNormalised(i,3) > ...
            xyzNormalised.imageStopNormalised(i,3)) 
        
        xyzNormalised.imageStartNormalised(i,3) = ...
            xyzNormalised.imageStartNormalised(i,3) + deltaxyz(i,3);
        
        xyzNormalised.imageStopNormalised(i,3) = ...
            xyzNormalised.imageStopNormalised(i,3) - deltaxyz(i,3);
    else
        xyzNormalised.imageStartNormalised(i,3) = ...
            xyzNormalised.imageStartNormalised(i,3) - deltaxyz(i,3);
        
        xyzNormalised.imageStopNormalised(i,3) = ...
            xyzNormalised.imageStopNormalised(i,3) + deltaxyz(i,3);
    end
    
end
    xyzLength = xyzLength + 2*deltaxyz;
     
    xyzMid = (xyzNormalised.imageStopNormalised + ...
    xyzNormalised.imageStartNormalised)./2;

%increase Lengths to adjust for extra half pixels at beginning and end
  
 end

