figure()
hold on;
L = 5e-2;
V = 613;
lambda = 800e-9;
colour = ['r';'b';'g']
colourNow = 0;
for F = [-1 1e9 1]  
    
    dwellTimeNs = [25,30,35,40,45,50,100,150,200,300,400,600,800,1000,1200,1400,1600];
    scanTimeNs = dwellTimeNs * 200;
    
    angVel = 17.2e6 ./ scanTimeNs;
    v = angVel * F;
    A = -v / V;
    
    firstChirp = V^2/lambda * (1 + A)./( (1 + A)*2*L + 2*(L + F) );
    secondChirp = V^2/lambda * (1 - A)./(2 * (L + F) );
    
    absDiff(1,:) = abs(firstChirp .* scanTimeNs * 1e-9);
    absDiff(2,:) = abs(secondChirp .* scanTimeNs * 1e-9);
    absDiff(3,:) = abs(secondChirp .* scanTimeNs * 1e-9);
    
    colourNow = colourNow + 1;
    
    plot(dwellTimeNs, absDiff(colourNow,:)*1e-6, colour(colourNow));
    ylim([0,80])
    xlabel('dwell time / ns')
    ylabel('start-stop absolute frequency difference / MHz')
end
legend([{' -100 \mum'} {'     0 \mum'} {'+100 \mum'}])
hold off;