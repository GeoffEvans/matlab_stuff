
%sawtooth test uses Hexfun 26th Nov 13
% Orbit_velocity adds in a correction for the velocity of the focal point
% see notes 11th Nov13
%..Orbit_position has only the position component of frequency deviation in
%the equation for frequency, it none the less produces a focussed beam.
%Varying the deflection fraction Fr does not change focus only magnitude of
%deflection
%oscilating along a line close to y=0 at the correct frequency, it is not in the correct place, i.e. it
%is at infinity as shown by the focus lens
%Orbit 5th November adds time varying frequency drives for the first time
%AOL6rays  is a program based on notes 14th Oct 13 aimed at modelling
%the ray paths of rays through a 6 AOD AOL
% Assume light propagating in in -Z direction and input plane is plane 1 so
% AODs are numbered 2:2+N-1 where N is the number of AODs
% Assume initially that input and AOD planes are planar and normal to Z at
% positions Zp. I is Z position of the image plane 
close all
clear 
%%%%%%%%%%%%%%%%%%%%%
deltaf= 1e6; %maximum frequency deviation
a= 10e-6; % period of hexfun
tc=[1 2 3 4 5 6 ].* a./6 ; % time sampling points, not necessarily evenly spaced, tc(6)=a
fin=[ -deltaf -deltaf 0 deltaf deltaf 0]; % frequency into sampling point at end of period(n)
fout=[ -deltaf -deltaf 0 deltaf deltaf 0]; % frequency out of sampling point at beginning of period(n+1)
desc=[tc;fin;fout];
%%%%%%%%%%%%%%%%%%%%%%%
lambda0=0.8e-6; % wavelength
% now set up angles of the 6 AODs
Vac=619;  %acoustic velocity m/s
n=1:6;  % AOD number
%f= 35e6;  %Frequency of each AOD
mode=-1; %diffraction order of the AODs
thetan=(4*(-n+2)+0)*pi./12; % XY plane angle of AOD wave vector wrt X axis, 
% RH coordsstarts at 60 degrees then clockwise in steps of 60 degrees
Ioz=0; % Z coordinate of input plane
        Zaod= [1 2 3 4 5 6 ]*-50e-6; % NOTE 50e-3 normally Z coordinates of AOD planes
Imz= [ -.4 -1.4 -1.5]; % Z coordinates of image plane(s)
Zp=[Ioz Zaod Imz]; % concatenated list of the Z coordinates of all the Z 
% planes where rays can start, change direction or be displayed
% make Aperture Input array use complex numbers to represent XY plane
Aperture=15e-3;  % Aperture diameter mm
T=  [-10:1:10]*1e-6;  % Time steps
Arr=7;  %number of radial points in aperture
x = -Arr:Arr;
[X,Y,time] = meshgrid(x,x,T);
Array= Aperture/2*(X.*exp(1i.*pi./3)+ Y.* exp(1i.*0));
Array(abs(Array)>Aperture/2*Arr*3^0.5/2*0.99)=NaN;
Rays=ones(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2));

Rays(1,:,:,:) = 1.2* Array./Arr; %scales position vectors of the rays in the first input plane
Deflect=zeros(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2)); 
%Forbit=0;  % 1.0e3; % freqency of the orbit of the focal point Hz
Fr=1/3; % Fraction of total deflection provided by each AOD
F=zeros(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2)); %sets up array of frequencies at every point where a ray intersects an AOD or image plane
Fc=40e6; % central frequency of AODs

RayDirVec= zeros(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2)); %Rays(1,:,:)/10; % Ray Direction vector angle in radians 
% for a paraxial lens of 1m focal length at input plane
% note this is also the differential displacement vector array that would occur at unit
% Z displacement (1m)again the complex number represents the X ,Y
% angles/displacements@unit Z 
% see notes 5th Nov 13

% THETA=pi/5*th; % angle of focal point from XY=0 in focal XY plane wrt to X (longitude)  (radians)
% PHI = -10e-3; % angle of focus from origin wrt Z axis Latitude up from south pole radians
% Zf= -1.4e0+((ph-1)*0.4);  %Z coordinate of intended focus
% Rxy=Zf.*tan(PHI);
% DirF=zeros(size(n,2),size(T,2));
% l=zeros(size(Rays));
% targetx=real(Zp(9)*PHI.*exp(2.*pi.*1i.*THETA));
% targety=imag(Zp(9)*PHI.*exp(2.*pi.*1i.*THETA));
%Rays(2,:,:,:)=Rays(1,:,:,:) may not be necessary?
for n=2:7
    Timearray(n,:,:,:)= time;
end

for n=2:7
    l(n,:,:,:)=(abs(Rays(1,:,:,:))).*cos(angle(Rays(1,:,:,:))-thetan(n-1));%  distance of each 
    %ray on each zplane from XY=0 resolved in direction of sound propagation
    % duplicated for each time
    Trays(n,:,:,:)=Timearray(n,:,:,:)-l(n,:,:,:)./Vac-(n-1).*a/3; % time reference for each ray
 for d1=1:size(Timearray,2)
     for d2=1:size(Timearray,3)
         for d3=1:size(Timearray,4)
    F(n,d1,d2,d3)=Hexfun(desc, squeeze(Trays(n,d1,d2,d3)));
         end
     end
 end
 Deflect(n,:,:,:)= mode.*F(n,:,:,:).*lambda0./Vac.*exp(1i*thetan(n-1)); % Ray Deflection vector caused by each Z plane of the AOL
end
%Calculate output ray position and direction from layers 2:7 = the AOL
for Zn=2:7   %size(Zp,2);
RayDirVec(Zn,:,:,:)= RayDirVec(Zn-1,:,:,:)+Deflect(Zn,:,:,:); % Ray direction vector justafter each Z plane
end
for Zn=1:7    %size(Zp,2)-1;
    Rays(Zn+1,:,:,:)=Rays(Zn,:,:,:)+RayDirVec(Zn,:,:,:).*(Zp(Zn+1)-Zp(Zn));
end
Deflect(8,:,:,:)= Rays(8,:,:,:).*0; % put in a lens at surface 8
% calculate rays from lens to end of system
for Zn=8:size(Zp,2);
RayDirVec(Zn,:,:,:)= RayDirVec(Zn-1,:,:,:)+Deflect(Zn,:,:,:); % Ray direction vector justafter each Z plane
end
for Zn=8:size(Zp,2)-1;
    Rays(Zn+1,:,:,:)=Rays(Zn,:,:,:)+RayDirVec(Zn,:,:,:).*(Zp(Zn+1)-Zp(Zn));
end


% for t=3
% for n=1:size(Zp,2)
% figure(n)
% plot(squeeze(Rays(n,:,:,t)),'.k')
% axis equal
% hold on
% end
% end
for n=1:size(Zp,2)
    Zarray(n,:,:)=ones(size(X,2),size(X,2)).*Zp(n);
end

for t=1:size(T,2);
figure
for n=1:size(X,2)
    for m=1:size(X,2)

plot3(squeeze(real(Rays(:,n,m,t))),squeeze(imag(Rays(:,n,m,t))),squeeze(Zarray(:,n,m)),'color','b')
axis([-0.012 0.012 -0.012 0.012 -1.6 0.1])
xlabel('x')
hold on
    end
end
end
  
break
figure(size(T,2)+1)
% % set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
% %       'DefaultAxesLineStyleOrder','-|--|:')
for nn=1:6

plot(T,squeeze(F(nn+1,8,8,:)),'color','b','LineWidth',nn)
hold on
end


figure(size(T,2)+2)
 plot(squeeze(real(Rays(9,:,8,:))),squeeze(imag(Rays(9,8,:,:))), 'x')
 title ('focal plane over all time')
 hold on
 plot(targetx,targety,'o')
 %XY plot of focal plane (plane 9) over time
 axis equal
 for t2= 10:14
 figure(size(T,2)+2+t2)
 plot(squeeze(real(Rays(9,:,:,t2))),squeeze(imag(Rays(9,:,:,t2))), 'x')
 end
t2 =15;
 figure(size(T,2)+3+t2)
 plot (squeeze(real(Rays(9,8,8,:))),squeeze(imag(Rays(9,8,8,:))), 'x');
 title ('central ray for each time at focal plane -1.4m')
  figure(size(T,2)+4+t2)
  axis square
  axis equal
 plot (squeeze(real(Rays(8,8,8,:))),squeeze(imag(Rays(8,8,8,:))), 'x');
  figure(size(T,2)+5+t2)
  title ('central ray for each time at lens plane -400mm')
  axis equal
 plot (squeeze(real(Rays(7,8,8,:))),squeeze(imag(Rays(7,8,8,:))), 'x');
 axis equal
 title ('central ray for each time at last AOD -300mm')
 