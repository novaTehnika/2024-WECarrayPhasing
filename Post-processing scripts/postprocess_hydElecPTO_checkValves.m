%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_hydElecPTO_checkValves")
        load(files(j).name,'-regexp','^(?!out)\w')

        eff_wecPump_array(iVar,SS) = eff_wecPump;
        p_min_wp_array(iVar,SS) = p_min_wp;

    end

end

%% Find indices for missing data files

files = dir;
nfiles = size(files,1);
Done = [];
notDone = 1:50*5;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_hydElecPTO_checkValves")
        load(files(j).name,'-regexp','^(?!out)\w')
        [r,c,val] = find(notDone==iVar);
        notDone = [notDone(1:c-1), notDone(c+1:end)];
        Done = [Done, iVar];

    end

end

try 
    doneArrayStr = num2str(Done(1));
    for j = 2:length(Done)
        doneArrayStr = append(arrayStr,[',',num2str(Done(j))]);
    end
catch
    % just move on
end

try
    jobArrayStr = num2str(notDone(1));
    for j = 2:length(notDone)
        jobArrayStr = append(jobArrayStr,[',',num2str(notDone(j))]);
    end
    
    
    if 1
    for j = 1:length(notDone)
        iVar = notDone(j);
        eff_wecPump_array(iVar) = nan;
        p_min_wp_array(iVar) = nan;

    end
    end

catch
    % just move on
end

%% Plot min Pressure as a function of total accumulator volume
SS = 2;

black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

linestyles = {'-', '--', ':', '-.'};

supTitleFontSize = 9;
subTitleFontSize = 9;
axFontSize = 8;
bottomEdge = 1;
leftEdge = 3;
width = 86/25.4; % one column: 3+9/16, two column: 7.5
height = 3;
lineWidth = 0.5;
fontName = 'Palatino';

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];
set(fig,'defaultAxesColorOrder',[black; black]);

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = fontName;
ax1.FontSize = axFontSize;


xlabel('flow coefficient, low-pressure (L/s/kPa\textsuperscript{1/2})', ...
'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)
title(['WEC-Driven Pump Performacne Vs.',newline, ...
    'Check Valve Sizing: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname',fontName)

yyaxis left
semilogx(X*kv*1000*sqrt(1000),eff_wecPump_array(:,SS),'k-')
hold on
ylabel('WEC-driven pump efficiency', ...
'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)
ylim([0.7 1])

yyaxis right
hold on
semilogx(X*kv*1000*sqrt(1000),1e-3*p_min_wp_array(:,SS),'k--')
ylabel('minimum pressure in pump (kPA)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)

grid on

% xlim([1 100])

ax = gca;
ax.FontName = fontName;
ax.FontSize = axFontSize;

leg = legend('pump efficiency','min. pressure');
leg.FontSize = axFontSize;
leg.FontName = fontName;
set(leg, 'Location', 'best')