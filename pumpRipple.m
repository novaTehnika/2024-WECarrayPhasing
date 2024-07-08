t = linspace(0,3,1e3)';

figure
hold on
for nPump = 1:10
    pumpTS(:,nPump) = abs(sin(2*pi*t));
    phaseShift_perPump = pi/nPump;
    for i = 2:nPump
        pumpTS(:,nPump) = pumpTS(:,nPump) + ...
            abs(sin(2*pi*t + (i-1)*phaseShift_perPump));
    end
    plot(t,pumpTS(:,nPump)/nPump)
    flowRipple(nPump) = range(pumpTS(:,nPump))/mean(pumpTS(:,nPump));
end

figure
scatter(1:numel(flowRipple),flowRipple)