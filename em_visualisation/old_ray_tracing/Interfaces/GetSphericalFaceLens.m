function lens = GetSphericalFaceLens(position, focalLength)
  global maxY;
  shiftPos = sqrt(focalLength^2 - maxY^2);
  l1 = Interface(position + shiftPos, @(Y) CircularConcave(focalLength,Y), 1.5);
  l2 = Interface(position - shiftPos, @(Y) CircularConvex(focalLength,Y), 1);
  lens = [l1 l2];
end


