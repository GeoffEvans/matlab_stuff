function [  ] = efficiency_analytic()

nExt = 2.25985;
nExtCoeff = 2*2.24e-5*(180/pi)^2;
nOrd = 2.25955;
nOrdCoeff = 2*1.96e-5*(180/pi)^2;
k = 800e-9;
aodDepth = 0;

    kDiffUnit
    kInUnit
    acFreqs = 
    waveVectorsUnit = [];
    waveVectorsR = waveVectorsUnit;
    normals = waveVectorsUnit;

    for n = 1:len(waveVectorsUnit)
        normals(n) = Rotation(angles(n)) * [0;0;1];
        waveVectorsR(n) = Rotation(angles(n)) * waveVector(n);
        K(n) = acFreqs(n) * waveVectorsR(n);
        momentumDiff = k*(RefIndOrd()*kDiffUnit - RefIndExt()*kInUnit) - K(n);
        momentumDiffNormal(n) = dot(normals(n),momentumDiff);
        transmittance(n) = sinc(momentumDiffNormal(n)*aodDepth/2/pi);
    end    
    
    efficiency = prod(transmittance)^2;
end

function r = Rotation(angles)
    angleX = angles(1);
    angleY = angles(2);
    r = eye(3);
    r(1,3) = angleX;
    r(3,1) = -angleX;
    r(2,3) = angleY;
    r(3,2) = -angleY;
end

function nd = RefIndOrd(cosine)
    nd = nOrd + nOrdCoeff*(1-cosine);
end

function nd = RefIndExt(cosine)
    nd = nExt + nExtCoeff*(1-cosine);
end