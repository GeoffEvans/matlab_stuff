function [ X2,Y2,X3L,Y3L,X4L,Y4L,X3R,Y3R,X4R,Y4R,lPrh,rPrh ] = CalculateCoords( xMax,yMax,O,S )
%Calculate the coords of the polarisations in the different plots
%Uses phase = wt +/- kx - w*timeLag

phaseDif = pi/2; % for principal ellipse

% Calculate rotated ellipse by rotation. This is for left handed ellipse.
x2 = xMax*cos(O) - yMax*exp(-1i*phaseDif)*sin(O);
y2 = xMax*sin(O) + yMax*exp(-1i*phaseDif)*cos(O);

% Find new ellipse params
x2Max = abs(x2);
y2Max = abs(y2);
phase2DifferenceL = angle(y2) - angle(x2);
phase2DifferenceR = -phase2DifferenceL; % arg(z') = -arg(z)
X2 = real(x2Max * exp(1i * S));
Y2 = real(y2Max * exp(1i * (S - phase2DifferenceL)));

% y is slow axis so phase of y will lag behind phase of x
% phaseY = wt +/- kx - phaseDif
% phaseY therefore incurs a negative penalty (lagging so needs longer to get the same phase)
waveplate = 1/4;
phase3DifferenceL = phase2DifferenceL - 2*pi*waveplate; 
phase3DifferenceR = phase2DifferenceR - 2*pi*waveplate;

X3L = real(x2Max * exp(1i * S));
Y3L = real(y2Max * exp(1i * (S - phase3DifferenceL)));
X3R = real(x2Max * exp(1i * S));
Y3R = real(y2Max * exp(1i * (S - phase3DifferenceR)));

X4L = X3L * cos(O) + Y3L * sin(O);
Y4L = Y3L * cos(O) - X3L * sin(O);
X4R = X3R * cos(O) + Y3R * sin(O);
Y4R = Y3R * cos(O) - X3R * sin(O);

% for phaseDifference of pi/2 we have right handed. 
lPrh = mod(phase3DifferenceL,2*pi) < pi;
rPrh = mod(phase3DifferenceR,2*pi) < pi;

end

