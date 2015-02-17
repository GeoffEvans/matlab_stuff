function TestAngleFunctionsAod3d()
TestMakeThetasPositive();
TestChangingAxis();
end

function TestMakeThetasPositive()
t=  rand(5) * 100 - 50;
p = rand(5) * 100 - 50;
[t1,p1]=aod3d.MakeAllThetaPositive(t,p);
shiftP = p1 - mod(p,2*pi);
negativeT = t < 0;
mod((negativeT-1).*shiftP, 2*pi) % Ans 1
v = get_vector_from_angles(1,t,p);
v1 = get_vector_from_angles(1,t1,p1);
vDiff = abs(v - v1) > 1e-13 % Ans 2
end

function TestOpticalRotationAngles()
thetaIn = 0:0.001*pi:2*pi;
opticAngle = pi/4;
[ theta, phi] = aod3d.OpticalRotationAngles(opticAngle, thetaIn);
vRotated = get_vector_from_angles(1,theta,phi);
plot3((vRotated(1,:)+vRotated(2,:))/sqrt(2), (-vRotated(1,:)+vRotated(2,:))/sqrt(2), vRotated(3,:));
hold on;
v = get_vector_from_angles(1,thetaIn,pi/4);
plot3((v(1,:)+v(2,:))/sqrt(2), (-v(1,:)+v(2,:))/sqrt(2), v(3,:), 'r');
vAngle = get_vector_from_angles(1,0:0.01*pi:opticAngle,pi/4+pi/2);
plot3((vAngle(1,:)+vAngle(2,:))/sqrt(2), (-vAngle(1,:)+vAngle(2,:))/sqrt(2), vAngle(3,:), 'g');
grid on;
axis square;
xlabel('<110>')
ylabel('<~10>')
zlabel('<001>')
hold off;
end

function TestChangingAxis() % Needs thetaOffset or phiOffset set to zero.
thetaIn = 0.9;
phiIn = 8;
thetaOffset = 0;
phiOffset = pi/2;
thetaX = thetaOffset+pi/2; % 30'
phiX = phiOffset;
thetaZ = thetaOffset;
phiZ = phiOffset;
[ theta, phi] = aod3d.ConvertAnglesBetweenFrames(thetaIn, phiIn, thetaZ, phiZ, thetaX, phiX)
thetaExpected = thetaIn - thetaOffset;
phiExpected = phiIn - phiOffset;
[thetaExpected,phiExpected] = aod3d.MakeAllThetaPositive(thetaExpected,phiExpected) 
end