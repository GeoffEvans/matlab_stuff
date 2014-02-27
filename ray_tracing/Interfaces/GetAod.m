function lens = GetAod(positionOfCentre, centreFreq, chirp, time)
% Take wavelength = Vac = 1
focalDepth = 1/chirp;
deflection = centreFreq + chirp * time;
l1 = GetWedge(positionOfCentre, deflection * 2);
l2 = GetPlanoConvex(positionOfCentre, focalDepth);
lens = [l1 l2];
end


