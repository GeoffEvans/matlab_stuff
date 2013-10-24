% For various wavevector directions, l, find the phase velocity eigenvalues 
% and displacement velocity eigenvectors. 

% Orientations of input plane waves
thetaRange = 0:pi/90:pi;
phiRange = 0:pi/90:2*pi;
[thetaMesh, phiMesh] = meshgrid(thetaRange, phiRange);
thetaPhiMesh = arrayfun(@(x,y) [x,y], thetaMesh, phiMesh, 'UniformOutput', false);

[ velocities, eVecs ] = find_v_for_angle(thetaMesh, phiMesh);

% Calculate polarisation from eigenvectors
HowTransverse = @(t,p,v,i) norm(cross([sin(t).*cos(p) sin(t).*sin(p) cos(t)], v(:,i)));
HowTransverse1 = @(t,p,v) HowTransverse(t,p,v,1);
HowTransverse2 = @(t,p,v) HowTransverse(t,p,v,2);
HowTransverse3 = @(t,p,v) HowTransverse(t,p,v,3);
c1 = cellfun(HowTransverse1, num2cell(thetaMesh), num2cell(phiMesh), eVecs);
c2 = cellfun(HowTransverse2, num2cell(thetaMesh), num2cell(phiMesh), eVecs);
c3 = cellfun(HowTransverse3, num2cell(thetaMesh), num2cell(phiMesh), eVecs);

% Get coords for plotting
GetK = @(i) cellfun(@(v) 1./v(i), velocities);
k1Mesh = GetK(1);
k2Mesh = GetK(2);
k3Mesh = GetK(3);

X1 = k1Mesh .* sin(thetaMesh).* cos(phiMesh);
Y1 = k1Mesh .* sin(thetaMesh).* sin(phiMesh);
Z1 = k1Mesh .* cos(thetaMesh);
X2 = k2Mesh .* sin(thetaMesh).* cos(phiMesh);
Y2 = k2Mesh .* sin(thetaMesh).* sin(phiMesh);
Z2 = k2Mesh .* cos(thetaMesh);
X3 = k3Mesh .* sin(thetaMesh).* cos(phiMesh);
Y3 = k3Mesh .* sin(thetaMesh).* sin(phiMesh);
Z3 = k3Mesh .* cos(thetaMesh);

s1 = surf(X1,Y1,Z1,c1,'EdgeColor','none');
hold on;
s2 = surf(X2,Y2,Z2,c2,'EdgeColor','none');
s3 = surf(X3,Y3,Z3,c3,'EdgeColor','none');
% set(gcf, 'renderer', 'zbuffer'); % Needed for colorbar to work on ucbtgje
% % but breaks alpha.
grid on;
axis square;
xlabel('x');
ylabel('y');
zlabel('z optic axis');
alpha(s1, 0.2);
alpha(s2, 0.2);
alpha(s3, 0.2);
alphamap('rampdown'); 
camlight right; 
lighting phong;
