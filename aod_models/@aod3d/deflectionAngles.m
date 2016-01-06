% All angles in radians
f = 20:0.01:60;
SimpleModelDeflection = @(f) 800e-9 .* f * 1e6 / 613 * 180/pi; % Simple model
clf
plot(f,SimpleModelDeflection(f));
xlabel('acoustic frequency / MHz')
ylabel('deflection angle / degrees')
hold on;
ai = 0.3425;
af = 25.6233;
at = 342.4658;
NewModelDeflection = @(f,i) -(1 - af * f - ai * i) / at; % data from Maak model
plot(f,NewModelDeflection(f,0.1), 'r');
PreliminaryResultsDeflection = @(f) (10.5/15 * (f-25) + 17)/527*180/pi;
plot(f,PreliminaryResultsDeflection(f), 'g');