trials = 10;
points = 15;
pointsArray = 1:points;
timeData = zeros(6,points,trials);
baseMultiple = 5000;
for M = 1:trials
    for N = pointsArray
        N %Print
        count = baseMultiple * N;
        [timeData(1,N,M),timeData(2,N,M),timeData(3,N,M),timeData(4,N,M),timeData(5,N,M),timeData(6,N,M)] = SpeedTest(count);
    end
end
% Plot
figure();
hold on;
lengthsArray = baseMultiple * pointsArray;
colourMap = hsv(7);
timeAverage = sum(timeData, 3) / trials;
for k = 1:6
    plot(lengthsArray,timeAverage(k,:),'color',colourMap(k,:));
end
legend('arrayfun','for-loop (array)','for-loop (cell)','100 * premultiplied','100 * MEX','parfor');
xlabel('runs')
ylabel('time / s')
hold off;