function [tRamp,A,B,C] = chirpCompute(scanVar,systemVar,xyzNormalised)

[ A, B, C, tRamp ] = driver_wrapper( scanVar,systemVar,xyzNormalised );

%  if (scanVar.aolMode == 1)
%     [tRamp, A, B, C] = miniScan( scanVar,systemVar,xyzNormalised);
%  end
 
% 
% [ A, B, C, tRamp ] = driver_wrapper( scanVar,systemVar,xyzNormalised );
% end
% 
%  if (scanVar.aolMode == 2)
% [ A, B, C, tRamp ] = driver_wrapper( scanVar,systemVar,xyzNormalised );
%  
%  end
% 
% % if aolMode == 3
% % 
% % [Tramp, A, B, C] = functionalImage( ~ )
% % 
% % if aolMode == 4
% % [Tramp, A, B, C] = miniScan( ~ )
end

