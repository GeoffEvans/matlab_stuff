range = 0:40;
modelVals = range;
nOrd = range;
nExt = range;
v_ac = range;
for j = range;
    angle = j/max(range)*5*pi/180;
    modelVals(j+1) = phase_matching(angle);
    [ nOrd(j+1), nExt(j+1), ~, ~, ~ ] = find_n_for_angle(angle);
    v_ac(j+1) = cellfun(@min, find_v_for_angle(pi/2-angle, pi/4));
end

degRange = range/max(range)*5;

plot(degRange, - modelVals * 180 / pi, 'color', 'red');
hold on;
% 
% acFreq = 35e6;
% opWavelenVac = 800e-9; 
% eqMaxVals = opWavelenVac * acFreq ./ v_ac ./ nOrd;
% eqMinVals = opWavelenVac * acFreq ./ v_ac ./ nExt;
% plot(degRange, eqMinVals * 180 / pi, 'color', 'blue');
% plot(degRange, eqMaxVals * 180 / pi, 'color', 'green'); 
xlabel('incidence / degrees');
ylabel('deflection / degrees');