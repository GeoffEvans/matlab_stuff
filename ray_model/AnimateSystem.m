% Animates lenses with "opposite" refractive indices
% Sweeps up, down, right, inf, left.

global opticalLength maxY minY stepY;
opticalLength = 200;
maxY = 10;
minY = -10;
stepY = 0.01;

fig = figure();
title('Lenses','Color',[.6 0 0]);

% Control smoothness of animation
k1 = 50;
k2 = 50;
k3 = 200;
k4 = 50;
k5 = 50;
k6 = 50;
count = 0;

for k = 1:k1
    rays{k} = PointSourceRays([0, 5*k/k1], 100, 1); % move up
end
count = count + k1;
for k = 1:k2
    rays{count+k} = PointSourceRays([0, 5 - 5*k/k2], 100, 1); % back down
end
count = count + k2;
for k = 1:k3
    rays{count+k} = PointSourceRays([200*k/k3, 0], 100, 1); % move to the right
end
count = count + k3;
for k = 1:k4
    rays{count+k} = PointSourceRays([200 + exp(exp(k4)), 0], 100, 1); % out to inf
end
count = count + k4;
for k = 1:k5
    rays{count+k} = PointSourceRays([-200-exp(50*(1-k/k5)), 0], 100, 1); % come round the left
end
count = count + k5;
for k = 1:k6
    rays{count+k} = PointSourceRays([-200*(1-k/k6), 0], 100, 1); % back to start
end
count = count + k6;

interfaces1 = GetSphericalFaceLens(100,40);
interfaces2 = GetInvertedSphericalFaceLens(100,40);

subplot(2,1,1);
ViewSystem(interfaces1, rays{1});
subplot(2,1,2)
ViewSystem(interfaces2, rays{1});
f(1) = getframe(gcf);
[im,map] = rgb2ind(f(1).cdata, 256, 'nodither');
 
for j = 2:count
    subplot(2,1,1, 'replace');
    ViewSystem(interfaces1, rays{j});
    subplot(2,1,2, 'replace')
    ViewSystem(interfaces2, rays{j});
    f(j) = getframe(gcf);
    im(:,:,1,j) = rgb2ind(f(j).cdata, map, 'nodither');
end

imwrite(im,map, strcat('animation', '.gif'), 'DelayTime', 0, 'LoopCount', inf)

    