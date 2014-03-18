
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
lambda0=0.8e-6; % wavelength
% now set up angles of the 6 AODs
Vac=619;  %acoustic velocity m/s
n=1:6;  % AOD number
%f= 35e6;  %Frequency of each AOD
mode=-1; %diffraction order of the AODs
thetan=(4*(-n+2)+0)*pi./12; % XY plane angle of AOD wave vector wrt X axis,
% RH coordsstarts at 60 degrees then clockwise in steps of 60 degrees
Ioz=0; % Z coordinate of input plane
Zaod= [1 2 3 4 5 6 ]*-50e-3; % Z coordinates of AOD planes
Imz= [ -.4 -0.6]; % Z coordinates of image plane(s)
Zp=[Ioz Zaod Imz]; % concatenated list of the Z coordinates of all the Z
% planes where rays can start, change direction or be displayed
% make Aperture Input array use complex numbers to represent XY plane
Aperture=15e-3;  % Aperture diameter mm
T = [0 1e-4]% [-10:1:10]*0.2e-5;  % Time steps
Arr=7;  %number of radial points in aperture
x = -Arr:Arr;
[X,Y,time] = meshgrid(x,x,T);
Array= Aperture/2*(X.*exp(1i.*pi./3)+ Y.* exp(1i.*0));
Array(abs(Array)>Aperture/2*Arr*3^0.5/2*0.99)=NaN;
Rays=ones(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2));

Rays(1,:,:,:) = 1.2* Array./Arr; %scales position vectors of the rays in the first input plane
DeflectOnAodN=zeros(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2));
Forbit=0;  % 1.0e3; % freqency of the orbit of the focal point Hz
Fr=1/3; % Fraction of total deflection provided by each AOD
F=zeros(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2)); %sets up array of frequencies at every point where a ray intersects an AOD or image plane
Fc=40e6; % central frequency of AODs

RayDirVecAfterAodN= zeros(size(Zp,2),2*Arr+1,2*Arr+1,size(T,2)); % Ray Direction vector angle in radians

THETA=pi/5; % angle of focal point from XY=0 in focal XY plane wrt to X (longitude)  (radians)
PHI = -10e-3; % angle of focus from origin wrt Z axis Latitude up from south pole radians
Zf= -0.5e0;  %Z coordinate of intended focus
Rxy=Zf.*tan(PHI);
DirF=zeros(size(n,2),size(T,2));
l=zeros(size(Rays));
targetx=real(Zp(9)*PHI.*exp(2.*pi.*1i.*THETA));
targety=imag(Zp(9)*PHI.*exp(2.*pi.*1i.*THETA));
%Rays(2,:,:,:)=Rays(1,:,:,:) may not be necessary?
for n=2:7
    Timearray(n,:,:,:) = time;
end

for n=2:7
    l(n,:,:,:)=(abs(Rays(1,:,:,:))).*cos(angle(Rays(1,:,:,:))-thetan(n-1));
    % distance of each ray on each zplane from XY=0 resolved in direction of sound propagation duplicated for each time
    Trays(n,:,:,:)= Timearray(n,:,:,:)-l(n,:,:,:)./Vac; % time reference for each ray
    F(n,:,:,:)=1./3.*Vac.*(((Zf-Zp(n)).^2+Rxy.^2).^-0.5).*Trays(n,:,:,:); % radial focus component of frequency
    DeflectOnAodN(n,:,:,:)= mode.*F(n,:,:,:).*exp(1i*thetan(n-1)); % Ray Deflection vector caused by each Z plane of the AOL
end
%Calculate output ray position and direction from layers 2:7 = the AOL
for Zn=2:7   %size(Zp,2);
    RayDirVecAfterAodN(Zn,:,:,:)= RayDirVecAfterAodN(Zn-1,:,:,:)+DeflectOnAodN(Zn,:,:,:); % Ray direction vector justafter each Z plane
end
for Zn=1:7    %size(Zp,2)-1;
    Rays(Zn+1,:,:,:)=Rays(Zn,:,:,:)+RayDirVecAfterAodN(Zn,:,:,:).*(Zp(Zn+1)-Zp(Zn));
end
for Zn=8:size(Zp,2)-1;
    Rays(Zn+1,:,:,:)=Rays(Zn,:,:,:)+RayDirVecAfterAodN(7,:,:,:).*(Zp(Zn+1)-Zp(Zn));
end

%%%%%%%% Plotting only below
figure()
hold on;
for Zn=1:9    %size(Zp,2)-1;
    fill3([Arr Arr -Arr -Arr],[Arr -Arr -Arr Arr],repmat([Zp(Zn)],1,4),repmat([Zp(Zn)],1,4))
end
alpha(0.1)
for n=1:size(Zp,2)
    Zarray(n,:,:)=ones(size(X,2),size(X,2)).*Zp(n);
end

for t=1:size(T,2);
    for n=1:size(X,2)
        for m=1:size(X,2)
            plot3(squeeze(real(Rays(:,n,m,t))),squeeze(imag(Rays(:,n,m,t))),squeeze(Zarray(:,n,m)),'color','b')
            axis([-0.012 0.012 -0.012 0.012 -0.6 0.1])
            xlabel('x')
            hold on
        end
    end
end
