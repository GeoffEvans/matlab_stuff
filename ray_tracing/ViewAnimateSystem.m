% Animates lenses with "opposite" refractive indices
% Sweeps up, down, right, inf, left.

global opticalLength maxY minY stepY;
opticalLength = 100;
maxY = 10;
minY = -10;
stepY = 0.01;

fig = figure();
title('Lenses','Color',[.6 0 0]);

chirp1 = 0.01;
chirp2 = 0.5 / (0.5/chirp1 - 20);
interfaces = [GetAod(20,0.04,chirp1,0) GetAod(40,-0.1,chirp2,0)];
rays = PlaneWaveRays(@(x) 1);
ViewSystem(interfaces, rays);
