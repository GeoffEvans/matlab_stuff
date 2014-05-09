
% 06_03_2014 - CHANGED to scale B (div by 8) as for scanning mode
% v4 corrects gridpoint for equations on patent 15th April 08
% v3 makes movex,movey and astig all functions of Z
% This program specifies the framework within which data is scanned, acquired and displayed.
% This is necessary so that aquired brightness data can be assigned unambiguously to a partiular voxel.
% In pointing mode this program makes a sequence of static points. for terminology and
% explanation of basics see Operating system for 3D v4.doc
% It is based on XYscan_atZ_seqMCPv2.m. Much of this is removed as this
% program does not scan, simply points for a specified time at particular
% points in 3D space. It is asssumed that a separate program will ensure
% that the points are within the accessible double pyramid shaped volume.
% A check ensures that RF drive does not go out of frequency range and
% notifies the user of any necessary reductions in dwell time or out of
% range points.
% it needs a variable pre loaded set of triggers to run the system
% All units are SI
%clear all
close all
format short e
lambda1=800;%%%%%%%%%%%H changed here
lambda = lambda1*1e-9;
%skewZX = 0
%skewZY = 0
Fc1 = 40; %%%%%%%%%%%%%H
Fc = Fc1*1e6;%%%%%%%%%H
Acc1 = 4.35; %%%%%%%%%%%%%H
Acc  = Acc1*1e-3;%%%%%%%%%H
sysclk = 240e6;%%%%%%%%%%%%H
daq_clock = 200e6;%%%%%%%%%%H
%CHANGED
%DATACLOCK

%dataclock_sys = (sysclk/25e6)/sysclk
%dataclock_daq = (daq_clock/25e6)/daq_clock
dataclock_sys = (sysclk/20e6)/sysclk; % This is modified to match 240MHz Clock rate
dataclock_daq = (daq_clock/20e6)/daq_clock; % This is modified to match 240MHz Clock rate


if (dataclock_sys == dataclock_daq)
    dataclock = dataclock_sys;    %  dwell time must be multipl of 50ns (removed *2)
else
    error('sysclock is not an absolute multiple of 25MHz')
end

% Make dwelltime an absolute multiple of dataclock
dwell = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%H
dwelltime = dwell*1e-6/dataclock
imaging =0;%%%%%%%%%%%%%%%H
skewZX=0;
skewZY=0;

zoom1 = 1;%%%%%%%%%H
zoom=round(zoom1) %  make 'zoom' an integer
Nvox1 = 100;%%%%%%%%%%%%%H
Nvox = round(Nvox1);
gridnumberXY= Nvox1; %%%%%%%%%%H
%
% if rem(Nvox,2)
%    gridnumberXY= Nvox1;
% else
%   gridnumberXY= Nvox1 - 1;
% end

gridnumberZ= 1;%%%%%%%%%%%%H
%lambda = 8e-7;

%Nvox=round(200) % integer number of unzoomed Voxels in each of the X,Y, and Z directions see diagram 10/03/07
% This is also the maximum dimension of any image
Mvox = Nvox*zoom; % integer number of zoomed voxels along one dimension

Xp=0;% X coordinat of centre of area to be imaged (Normalised to lie between -1 and +1)
Yp=0;% Y coordinate of centre of area to be imaged (Normalised to lie between -1 and +1)
Zp=0;
% the next set of 9 parameters fine tune the ramp rates of the 4 AODs to
% correct for the effect on position of the angled surfaces of he AODS.
% They are each of the form 1+alpha +beta*Znp + gamma* Znp^2
alphaX=0
betaX=0 % modifies MoveX
gammaX=0
alphaY=0
betaY=0 % modifies moveY
gammaY=0
alphaAstig=0
betaAstig=0  % modifies Astig, the astigmatism
gammaAstig=0

% Y coordinate of centre of area to be imaged (Normalised to lie between -1 and +1)
% an area will be scanned between Xc-1/zoom and Xc+1/zoom and Yc-1/zoom and
% Yc + 1/zoom. If some of the area comes outside the drive limits it will not
% be scanned, and on the final image its intensity will be put at zero.
%dataclock= 4e-6;
%triggerwidth=ceil(4e-6/dataclock)*dataclock;
%Amplitude  = ([2100 1500 1900 1500]); %Max acoustic drive amplitude for AODs 1-4, each in range 0-4095
% Final set drive Amplitude for a particular miniscan = Amplitude* AmpMultgp(k)(becomes AmplitudeM(p)), the Amplitude multiplier.
mode=-1 %diffraction order used by the AODs, changes sign of angle of deflection and wavefront curvature
%lambda=800e-9;
%NEW
%L=22e-3; % aperture of AOD m
L=15e-3; % aperture of AOD m
Va=619 ;% speed of sound m/s
AODfill=dataclock*floor(L/Va/dataclock); % the AOD fill time rounded up to an integral number of dataclock cycles

Deltafmax=Acc*Va/lambda ; %maximum permitted frequency deviation from Fc   Hz
d1=50e-3 ;%physical separation of plane AOD1(X1) from AOD2(Y1) m
AODth1 = 8e-3; % Thickness of TeO2 between AOD grating planes 1 & 2
d1app = (d1-AODth1+AODth1/2.3);
% apparent optical thickness between grating planes
d2=50e-3; %physical separation of plane AOD2 (Y1) from AOD3(X2) m 45micron opt 300mm
AODth2 =8e-3; % Thickness of TeO2 between AOD grating planes 2 & 3
d2app= (d2-AODth2+AODth2/2.3);
% apparent optical thickness between grating planes 2&3
d3=50e-3; %physical separation of plane AOD3 (X2) from AOD4 (Y2) m
AODth3 = 8e-3; % Thickness of TeO2 between AOD grating planes 3 & 4
d3app= (d3-AODth3+AODth3/2.3);
% apparent optical thickness between grating planes 3&4
d4appmin= Va^2*((L/Va))/(4*lambda*Deltafmax); %corrected 150408
% the minimum distance from AOD4 for to the apparent focus of the laser beam, It is
% positive and becomes smaller the greater the maximum possible ramp rate
%DATACLOCK
%dwelltime = floor (dwell)
% START OF POINTING MODE
% This first lists a set of Grid Point Centres(GPCs),then forms a 3D grid
% of points around each GPC
% Note in this program the subscripts p and m all refer to the same
% miniscan number. This for historical reasons
% for each of the points below will make a grid of points centred on it
Ngxy=round(gridnumberXY); % number of grid points per GPC in X and Y (must be odd)
Ngz=round(gridnumberZ);% number of grid points per GPC in Z(must be odd)
Spxy=round(1) %spacing of gridpoints in X and Y (voxels)
Spz=round(1) %spacing of gridpoints in Z(voxels)

AmpMultgp=ones(1,length(Xp))% amplitude multiplier for RF drive of each GP (for varying focussed spot intensity
% with e.g. depth(must be <1 or will cause overdrive of AODs). This
% multiplier is applied to all 4 AODs. As they are in series the intensity
% reduces approximatly as the fourth power of AmpMultgp
%CHANGED
%DATACLOCK
%NEW
Tdgp=dataclock*(dwelltime.*ones(1,length(Xp))); % dwell times for each point of the GPC sequence(=integer*dataclock  s )
%Tdgp=dataclock*ceil(dwelltime.*ones(1,length(Xp))); % dwell times for each point of the GPC sequence(=integer*dataclock  s )
%Tdgp=1e-6*ceil(dwelltime.*ones(1,length(Xp))); % dwell times for each point of the GPC sequence(=integer*dataclock  s )
NGPC=size(Zp,2) % number of Grid Point Centres per pointing mode cycle
NC=round(1) % number of pointing mode cycles per AOD for loading onto iDDS
if max(AmpMultgp)>1;
    error;('line56 AmpMultgp greater than 1')
end
PMCtime=dataclock*ceil(1e-3/dataclock) % Pointing mode cycle time (= integer*dataclock s )
Scantime=NC*PMCtime
% choose this to be slightly greater than sum of scan times of PMC,its
% function is to ensure monitored time points are at a convenient time
% interval. If sum of scan times  > PMCtime then PMC time becomes an increasing integer *
% PMCtime until PMCtime > sum of scan times.

% now adjust XYZ coordinates point precisely at cente of nearest voxel and
% save voxel XYZ identifiers as as a function of voxel number
for u=1:NGPC;  % first u loop finds voxel identity of each Grid Point Centre
    if (rem(Mvox,2))&(rem(Ngxy,2)); % if Mvox odd this statement is true (non zero)
        Mxgp(u)=floor(Xp(u)*(Mvox-1)/2)+Mvox/2 +0.5; %X zoomed voxel identifier (Mvox odd)
        Mygp(u)=floor(Yp(u)*(Mvox-1)/2)+Mvox/2 +0.5; %Y zoomed voxel identifier (Mvox odd)
        Mzgp(u)=floor(Zp(u)*(Mvox-1)/2)+Mvox/2 +0.5; %Z zoomed voxel identifier (Mvox odd)
        %changed for even number of pixels
        %elseif (rem(Ngxy,2))&(rem((Mvox+1),2));
    else
        %     Mxgp(u)=floor(Xp(u)*(Mvox/2))+Mvox/2+1;%X zoomed voxel identifier (Mvox even)
        %     Mygp(u)=floor(Yp(u)*(Mvox/2))+Mvox/2+1;%Y zoomed voxel identifier (Mvox even)
        %     Mzgp(u)=floor(Zp(u)*(Mvox/2))+Mvox/2+1;%Z zoomed voxel identifier (Mvox even)
        Mxgp(u)=floor(Xp(u)*(Mvox/2))+Mvox/2;%X zoomed voxel identifier (Mvox even)
        Mygp(u)=floor(Yp(u)*(Mvox/2))+Mvox/2;%Y zoomed voxel identifier (Mvox even)
        Mzgp(u)=floor(Zp(u)*(Mvox/2))+Mvox/2;%Z zoomed voxel identifier (Mvox even)
    end
end % end of first u loop
k=0

if imaging == 0
    for u=1:NGPC; % GPC loop
        for z=1:Ngz; % z grid loop
            for w=1:Ngxy; % y grid loop
                for v=1:Ngxy;% x grid loop
                    k=k+1;
                    Td(k)=Tdgp(u);
                    AmpMult(k)=AmpMultgp(u);
                    %changed for even number of pixels ???
                    %                 Mxp(k)=(Mxgp(u)+Spxy*(v-(Ngxy)/2));
                    %                 Myp(k)=(Mygp(u)+Spxy*(w-(Ngxy)/2));
                    %                 Mzp(k)=(Mzgp(u)+Spz*(z-(Ngz)/2));
                    
                    if rem(Mvox,2); % if Mvox odd this statement is true (non zero)
                        Mxp(k)=(Mxgp(u)+Spxy*(v-(Ngxy+1)/2));
                        Myp(k)=(Mygp(u)+Spxy*(w-(Ngxy+1)/2));
                        Mzp(k)=(Mzgp(u)+Spz*(z-(Ngz+1)/2));
                        Xnp(k)=2*(Mxp(k)-Mvox/2-0.5)/(Mvox-1);  % exact X coordinate of centre of pointing voxel
                        Ynp(k)=2*(Myp(k)-Mvox/2-0.5)/(Mvox-1);  % exact Y coordinate of centre of pointing voxel
                        Znp(k)=2*(Mzp(k)-Mvox/2-0.5)/(Mvox-1);  % exact Z coordinate of centre of pointing voxel
                        
                    else
                        Mxp(k)=(Mxgp(u)+Spxy*(v-(Ngxy)/2));
                        Myp(k)=(Mygp(u)+Spxy*(w-(Ngxy)/2));
                        Mzp(k)=(Mzgp(u)+Spz*(z-(Ngz)/2));
                        Xnp(k)=2*(Mxp(k)-Mvox/2-0.5)/(Mvox-1);  % exact X coordinate of centre of pointing voxel
                        Ynp(k)=2*(Myp(k)-Mvox/2-0.5)/(Mvox-1);  % exact Y coordinate of centre of pointing voxel
                        Znp(k)=2*(Mzp(k)-Mvox/2-0.5)/(Mvox-1);  % exact Z coordinate of centre of pointing voxel
                    end
                end
            end
        end
    end
    
else
    
    
    for u=1:NGPC; % GPC loop
        for z=1:Ngz; % z grid loop
            for v=1:Ngxy; % x grid loop
                k=k+1;
                %changed for even number of pixels ???
                %              Mxp(k)=(Mxgp(u)+Spxy*(v-(Ngxy)/2));
                %              Mzp(k)=(Mzgp(u)+Spz*(z-(Ngz)/2));  %NB this puts the z voxel at 50.5 which gives Znp=0
                %              Mxp(k)=(Mxgp(u)+Spxy*(v-(Ngxy+1)/2));
                %              Mzp(k)=(Mzgp(u)+Spz*(z-(Ngz+1)/2));
                Td(k)=Tdgp(u);
                AmpMult(k)=AmpMultgp(u);
                if rem(Mvox,2); % if Mvox odd this statement is true (non zero)
                    Mxp(k)=(Mxgp(u)+Spxy*(v-(Ngxy+1)/2));
                    Mzp(k)=(Mzgp(u)+Spz*(z-(Ngz+1)/2));
                    Xnp(k)=2*(Mxp(k)-Mvox/2-0.5)/(Mvox-1);  % exact X coordinate of centre of pointing voxel
                    Znp(k)=2*(Mzp(k)-Mvox/2-0.5)/(Mvox-1);  % exact Z coordinate of centre of pointing voxel
                else
                    %changed for even number of pixels ???
                    Mxp(k)=(Mxgp(u)+Spxy*(v-(Ngxy)/2));
                    Mzp(k)=(Mzgp(u)+Spz*(z-(Ngz)/2));      %NB this puts the z voxel at 50.5 which gives Znp=0
                    Xnp(k)=2*(Mxp(k)-Mvox/2-0.5)/(Mvox-1);
                    %Xnp(k)=2*(Mxp(k)-Mvox/2-0.5)/(Mvox);  % exact X coordinate of centre of pointing voxel
                    Znp(k)=2*(Mzp(k)-Mvox/2- 0.5)/(Mvox-1);  % exact Z coordinate of centre of pointing voxel
                end
            end
            k = 0
            for w=1:Ngxy; % y grid loop
                k = k+1;
                %changed for even number of pixels ???
                %Myp(k)=(Mygp(u)+Spxy*(w-(Ngxy+1)/2));
                
                if rem(Mvox,2); % if Mvox odd this statement is true (non zero)
                    Myp(k)=(Mygp(u)+Spxy*(w-(Ngxy+1)/2));
                    Ynp(k)=2*(Myp(k)-Mvox/2-0.5)/(Mvox-1);  % exact Y coordinate of centre of pointing voxel
                else
                    %changed for even number of pixels ???
                    Myp(k)=(Mygp(u)+Spxy*(w-(Ngxy)/2));
                    Ynp(k)=2*(Myp(k)-Mvox/2-0.5)/(Mvox-1);  % exact Y coordinate of centre of pointing voxel
                end
            end
        end
    end
    
end

NDV=size(Znp,2);  %Total number of Dwell voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next section computes the scan times and RF start and stop freqencies for
% each miniscan that points at a voxel of the sequence.
for k=1:NDV; % for each dwell voxel
    if Znp(k)==0;
        Znp(k)=0.0001;   %(to prevent divide by zero)
    end
    d4app=d4appmin/(Znp(k)); % this is the distance of the apparent focus from the final AOD.
    %It is positive for converging rays (Zn positive) and negative for
    %diverging rays
    %Next set the fine tuning of the ramp rate paramenters which are a function
    %of Znp
    movex=1+alphaX +betaX*Znp(k)+gammaX*Znp(k)^2
    movey=1+alphaY +betaY*Znp(k)+gammaY*Znp(k)^2
    astig =1+alphaAstig +betaAstig*Znp(k)+gammaAstig*Znp(k)^2
    % Next Calculate the four AOD ramp rates during pointing.
    a(1)=mode*((Va^2/lambda)/(2*(d4app+d3app)+d1app+d2app)...
        )/movex*astig;  %AOD X1 corrected 150408
    a(2)=mode*((Va^2/lambda)/(2*d4app+d2app+d3app)...
        )/movey/astig;   % AOD Y1 corrected 150408
    a(3)= mode*(Va^2/(2*lambda*(d4app+d3app)))* movex*astig; % AOD X2 corrected 150408
    a(4)=mode*(Va^2/(2*lambda*d4app))*movey/astig; % AODY2   corrected 150408
    
    %now calculate the X & Y frequency deflections to get the beam focussed on the
    %required point. See PK notes 27/07/07. Note all deflections are defined
    %wrt the centre of AOD4. First AOD on each axis deflects -ive frequency,
    %second, +ive. The frequency deflections scaled so that the anglular deflections are equal
    % and opposite at any fixed position on the second AOD of each pair.
    
    fdeflect(1)=-mode*2*(d3app+d4app)*d4app/((d4app+d3app)*(d1app+d2app+2*(d3app+d4app)))*(Xnp(k)+skewZX.*Znp(k))*Deltafmax;%minus
    % AOD X1 correct I think 150408
    fdeflect(2)=mode*2*d4app/(2*d4app+d3app+d2app)*(Ynp(k)+skewZY.*Znp(k))*Deltafmax;% AOD Y1 correct I think 150408
    fdeflect(3)=mode* d4app/(d4app+d3app)*(Xnp(k)+skewZX.*Znp(k))*Deltafmax;% AOD X2 correct I think 150408%minus
    fdeflect(4)=-mode* Ynp(k)*Deltafmax;% AOD Y2 correct I think 150408
    
    %NEW
    %AOD1 = X1 AOD2=Y1 AOD3 = X2 AOD4 = Y2 -- change the index for record
    % structure X1, X2, Y1, Y2
    timemid = -AODfill/2;
    
    fdeflect0(1) = fdeflect(1) + a(1)*timemid;
    fdeflect0(2) = fdeflect(2) + a(2)*timemid;
    fdeflect0(3) = fdeflect(3) + a(3)*timemid;
    fdeflect0(4) = fdeflect(4) + a(4)*timemid;
    
    
    A(k,1) = round((Fc + fdeflect0(1))*(2^32/sysclk));
    A(k,3) = round((Fc + fdeflect0(2))*(2^32/sysclk));
    A(k,2) = round((Fc + fdeflect0(3))*(2^32/sysclk));
    A(k,4) = round ((Fc + fdeflect0(4))*(2^32/sysclk));
    % Changed to scale like the raster for steep slopes
    Bfull(k,1) = round(a(1) * (2^32/sysclk^2));
    Bfull(k,3) = round(a(2) * (2^32/sysclk^2));
    Bfull(k,2) = round(a(3) * (2^32/sysclk^2));
    Bfull(k,4) = round(a(4) * (2^32/sysclk^2));
    %CHANGED to scale B as for scanning mode
    B(k,1) = round(Bfull(k,1)/8); % CHANGED JAN 2014
    B(k,2) = round(Bfull(k,2)/8);
    B(k,3) = round(Bfull(k,3)/8);
    B(k,4) = round(Bfull(k,4)/8);
    
    C(k,1) = 0;
    C(k,3) = 0;
    C(k,2) = 0;
    C(k,4) = 0;
    
    
    
    
    % Now calculate the maximum possible scan times for these combinations
    % of ramp rate and deflection
    %DATACLOCK- removed code for frequency limits
    %  DiffX=abs(fdeflect(1)-fdeflect(3));
    % tmaxX=2*(2*Deltafmax-DiffX)/(abs(a(1))+abs(a(3)));
    %DiffY=abs(fdeflect(2)-fdeflect(4));
    % tmaxY=2*(2*Deltafmax-DiffY)/(abs(a(2))+abs(a(4)));
    %if tmaxX<AODfill+Td(k);
    %  Td(k)=dataclock*floor(abs(tmaxX-AODfill)/dataclock);
    % Tdchd(k)=Td(k); % saves all changed dwell times (can even be negative
    % in which case the scan time will be less than the AOD fill time and
    % there will be no miniscan nor data acqisition
    %end
    %if tmaxY<AODfill+Td(k);
    %  Td(k)=dataclock*floor(abs(tmaxY-AODfill)/dataclock);
    %Tdchd(k)=Td(k); % saves all changed dwell times (can even be negative
    % in which case the scan time will be less than the AOD fill time and
    % there will be no miniscan nor data acqisition
    %end
    
    
    tscan(k) =AODfill+Td(k); % the scan time for each point
    
    
    Tramp(k,1)=  round(tscan(k)*sysclk);
    Tramp(k,2)=  round(tscan(k)*sysclk);
    Tramp(k,3)=  round(tscan(k)*sysclk);
    Tramp(k,4)=  round(tscan(k)*sysclk);
    
    fstart(k,1)=Fc-a(1)*tscan(k)/2+fdeflect(1);  % start frequency for AOD X1
    fstop(k,1)=Fc+a(1)*tscan(k)/2+fdeflect(1);   % stop frequency for AOD X1
    fstart(k,2)=Fc+fdeflect(2)-a(2)*tscan(k)/2;  % start frequency for AOD Y1
    fstop(k,2)=Fc+fdeflect(2)+a(2)*tscan(k)/2;    % stop frequency for AOD Y1
    fstart(k,3)=Fc-a(3)*tscan(k)/2+fdeflect(3);   % start frequency for AOD X2
    fstop(k,3)=Fc+a(3)*tscan(k)/2+fdeflect(3);    % stop frequency for AOD X2
    fstart(k,4)=Fc+0e6+fdeflect(4)-a(4)*tscan(k)/2;   % start frequency for AOD Y2
    fstop(k,4)=Fc+0e6+fdeflect(4)+a(4)*tscan(k)/2;    % stop frequency for AOD Y2
    
end %end of k loop
Totaldwellnumber=sum(Td)/dataclock*NC;
