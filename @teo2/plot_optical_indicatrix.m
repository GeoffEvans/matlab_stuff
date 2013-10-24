% For various wavevector directions, find the refractive indices and their respective polarisations

% Orientations of input plane waves
thetaRange = -pi/2:pi/360:pi/2;
% no phi dependence since paratellurite is uniaxial
[ n_ord, n_ext, p_ord, p_ext ] = find_n_for_angle(thetaRange);
activityVector = 2.65e-05; %See Warner White Bonner, 87deg/mm
f = figure();
subplot(1, 2, 1);
hold on;
grid on;
grid minor;
xlabel('theta');
ylabel('n');
thetaDegrees = thetaRange*180/pi;
plot(thetaDegrees,real(n_ord), 'color', 'red');
plot(thetaDegrees,real(n_ext));
delta = n_ord(1).^2 * activityVector / 2; %(1.83 Xu & St)
ordNumerator = ( cos(thetaRange)/((1 - delta)*n_ord(1)) ).^2 + ( sin(thetaRange)/n_ord(1) ).^2;
extNumerator = ( cos(thetaRange)/((1 + delta)*n_ord(1)) ).^2 + ( sin(thetaRange)/n_ext(1) ).^2;
nOrdApprox = ordNumerator .^ -0.5;
nExtApprox = extNumerator .^ -0.5;
plot(thetaDegrees,nOrdApprox, 'color', 'green');
plot(thetaDegrees,nExtApprox, 'color', 'magenta');
hold off;

% Calculate the handedness of the polarisations: negative phase implied
% left-handed, very big or small ratio implies linear
RightHanded = @(pol) angle(pol) > 1e-6;
Colour = @(pol) (RightHanded(pol)-0.5) .* log(abs(pol));

% Plot the possible n~ (will have two eigenvalues) and colour according to
% handedness.
phiRange = 0:pi/360:2*pi;
phiLength = length(phiRange);
[thetaMesh, phiMesh] = meshgrid(thetaRange, phiRange);
n1Mesh = ones(phiLength,1) * n_ord;
n2Mesh = ones(phiLength,1) * n_ext;
c1Mesh = ones(phiLength,1) * Colour(p_ord);
c2Mesh = ones(phiLength,1) * Colour(p_ext);

X1 = n1Mesh .* sin(thetaMesh).* cos(phiMesh);
Y1 = n1Mesh .* sin(thetaMesh).* sin(phiMesh);
Z1 = n1Mesh .* cos(thetaMesh);
X2 = n2Mesh .* sin(thetaMesh).* cos(phiMesh);
Y2 = n2Mesh .* sin(thetaMesh).* sin(phiMesh);
Z2 = n2Mesh .* cos(thetaMesh);

subplot(1, 2, 2);
s = surf(X1,Y1,Z1,c1Mesh,'EdgeColor','none');
% set(gcf, 'renderer', 'zbuffer'); % Needed for colorbar to work on ucbtgje
% % but breaks alpha.
alpha(s, 0.3);
grid on;
axis square;
hold on;
xlabel('x');
ylabel('y');
zlabel('z optic axis');
t = surf(X2,Y2,Z2,c2Mesh,'EdgeColor','none');
alpha(t, 0.2);

