function lens = GetInvertedSphericalFaceLens(position, focalLength)
  l1 = Interface(position - focalLength, @(Y) CircularConvex(focalLength,Y), 1.5);
  l2 = Interface(position + focalLength, @(Y) CircularConcave(focalLength,Y), 1);
  lens = [l1 l2];
end


