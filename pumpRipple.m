%% rotary piston pump flow ripple

t = linspace(0,3,1e3)';

figure
hold on
for nPump = 1:10
    pumpTS(:,nPump) = max(0,sin(2*pi*t));
    phaseShift_perPump = 2*pi/nPump;
    for i = 2:nPump
        pumpTS(:,nPump) = pumpTS(:,nPump) + ...
            max(0,sin(2*pi*t + (i-1)*phaseShift_perPump));
    end
    plot(t,pumpTS(:,nPump)/mean(pumpTS(:,nPump)))
    PDflowRipple(nPump) = range(pumpTS(:,nPump))/mean(pumpTS(:,nPump));
end

figure
scatter(1:numel(PDflowRipple),PDflowRipple)


%% WEC-driven pump flow ripple

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
    WPflowRipple(nPump) = range(pumpTS(:,nPump))/mean(pumpTS(:,nPump));
end

figure
scatter(1:numel(WPflowRipple),WPflowRipple)


%% Compared

f1 = figure;
scatter(1:numel(PDflowRipple),PDflowRipple,100,"black","x","LineWidth",1.5)
hold on
scatter(1:numel(WPflowRipple),WPflowRipple,100,"red","o","LineWidth",1.5)
xLim = xlim();
yLim = ylim();
clear f1

figure
semilogy([-1 -1],[1 1])
hold on
scatter(1:numel(PDflowRipple),PDflowRipple,100,"black","o","LineWidth",1.5)
scatter(1:numel(WPflowRipple),WPflowRipple,100,"red","x","LineWidth",1.5)

xlim(xLim);
% ylim(yLim);