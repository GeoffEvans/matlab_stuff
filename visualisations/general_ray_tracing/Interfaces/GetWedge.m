function lens = GetWedge(positionOfBack, angleFront)
  global maxY;
  shiftPos = abs(Tilt(maxY,angleFront));
  l1 = Interface(positionOfBack - shiftPos, @(Y) Tilt(Y,angleFront), 1.5);
  l2 = Interface(positionOfBack, @(Y) Y*0, 1);
  lens = [l1 l2];
end

function tilt = Tilt(Y, angle)
    tilt = Y .* tan(angle);
end


