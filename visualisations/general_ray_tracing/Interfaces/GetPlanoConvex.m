function lens = GetPlanoConvex(frontPosition, focalLength)
  global maxY;
  radius = focalLength/2;
  shiftPos = radius-sqrt(radius^2 - maxY^2);
  l1 = Interface(frontPosition + radius, @(Y) CircularConcave(radius,Y), 1.5);
  l2 = Interface(frontPosition + shiftPos, @(Y) Y*0, 1);
  lens = [l1 l2];
end


