%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_hydElecPTO_PnomHPaccum")
        load(files(j).name,'-regexp','^(?!out)\w')

        PP_WEC_array(iVar,SS) = PP_WEC;
        PP_wp_array(iVar,SS) = PP_wp;
        PP_mLoss_array(iVar,SS) = PP_mLoss;
        PP_gen_array(iVar,SS) = PP_gen;
        PP_hPRV_array(iVar,SS) = PP_hPRV;
        dpdt_max_array(iVar,SS) = dpdt_max;

    end

end

%% Find indices for missing data files

files = dir;
nfiles = size(files,1);
Done = [];
notDone = 1:nVar;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_hydElecPTO_PnomHPaccum")
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
        PP_WEC_array(iVar,SS) = NaN;
        PP_wp_array(iVar,SS) = NaN;
        PP_mLoss_array(iVar,SS) = NaN;
        PP_gen_array(iVar,SS) = NaN;
        PP_hPRV_array(iVar,SS) = NaN;
        dpdt_max_array(iVar,SS) = NaN;

    end
    end

catch
    % just move on
end


%% sort results into multi-dim. array

for iVar = 1:nVar
    iVar1 = find(p_nom == p_nom_mesh(iVar));
    iVar2 = find(V == V_mesh(iVar));
    PP_WEC(iVar1,iVar2,SS) = PP_WEC_array(iVar,SS);
    PP_wp(iVar1,iVar2,SS) = PP_wp_array(iVar,SS);
    PP_mLoss(iVar1,iVar2,SS) = PP_mLoss_array(iVar,SS);
    PP_gen(iVar1,iVar2,SS) = PP_gen_array(iVar,SS);
    PP_hPRV(iVar1,iVar2,SS) = PP_hPRV_array(iVar,SS);
    dpdt_max(iVar1,iVar2,SS) = dpdt_max_array(iVar,SS);
end


%% Plot power generated as a function of total accumulator volume
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





I = 1:2:numel(V);
for i = 1:numel(I)
    scatter(p_nom*1e-6,PP_gen(:,I(i),SS)*1e-3,100,color(i,:),'Marker','x','LineWidth',2)
    hold on
    legStr(i) = {[num2str(V(I(i))),' m^3']};
end
ylabel('Elec. power, mean (kW)', ...
    'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)
xlabel('Pressure (MPa)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)


grid on

% xlim([1 100])
title(['Mean Power Production Vs.', ...
    'Accumulator Volume:',newline, ...
    'Sea State ',num2str(SS),', ',...
    'Motor Disp. ',num2str(par.motor.D/1e-6*(2*pi)),'cc/rev'],...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname',fontName)

ax = gca;
ax.FontName = fontName;
ax.FontSize = axFontSize;

leg = legend(legStr);
leg.FontSize = axFontSize;
leg.FontName = fontName;
set(leg, 'Location', 'best')