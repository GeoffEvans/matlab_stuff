function [ inputRays ] = PlaneWaveRays( densityFunction )

% Need to recursively generate an array of differences
X = zeros(1,199);
Y = X;
angles = X;

for (i=1: length(Y)/2)
    Y(100+i) = Y(100 + i - 1) + 1 ./ densityFunction(Y(100 + i - 1));
    Y(100-i) = Y(100 - i + 1) - 1 ./ densityFunction(Y(100 - i + 1));
end

starts = ([X;Y]);
inputRays = Rays(starts, angles, []);

end

