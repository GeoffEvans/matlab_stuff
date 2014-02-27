function plotFourierSeries()
    x = -10:0.02:10;
    n = -50:50;
    fn = f(n);
    gn = g(n);
    hn = h(n);
    subplot(4,1,1);
    plot(x,fs(fn,n,x));
    axis([-10,10,-2,2])
    subplot(4,1,2);
    plot(x,fs(gn,n,x));
    axis([-10,10,-2,2])
    subplot(4,1,3);
    plot(x,fs(hn,n,x));
    axis([-10,10,-2,2])
    subplot(4,1,4);
    plot(x,fs(gn,n,x).*fs(fn,n,x));
    axis([-10,10,-2,2])
end

function fn = f(n)
    fn = -2*1i/pi./n;
    fn = fn .* mod(n,2); % even terms are zero
    fn(n==0) = 0; % no offset
end

function gn = g(n)
    gn = 1i/pi./n.*(-1).^n; 
    gn(n==0) = 0; % no offset
end

function hn = h(n)
    hn = -2/(pi.^2)./(n.^2);
    hn = hn .* mod(n,2); % even terms are zero
    hn(n==0) = 0.5; % half up
end

function ts = fs(cs,n,x) % assumes L = 2
    [mx,mc] = meshgrid(x,cs);
    [~,mn] = meshgrid(x,n);
    ts = mc .* exp(1i * pi .* mn .* mx);
    ts = sum(ts,1);
    ts = real(ts);
end
