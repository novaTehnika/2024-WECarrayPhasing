% To find PP_max, run the following from repo. root:
% load(["Sea States",filesep,"SSdata_Bull2017WEPrize.mat"])
% SS = 2;
% PP_max = max(PP_w_data(SS,:));
PP_max = 3.414e+05;

%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_hydElecPTO_arrayHPaccum")
        load(files(j).name,'-regexp','^(?!out)\w')

        PP_WEC_array(iVar,SS) = sum(PP_WEC);
        PP_wp_array(iVar,SS) = sum(PP_wp);
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
    if strfind(files(j).name,"data_hydElecPTO_arrayHPaccum")
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


%% sort results into multi-dim. array with optimal nominal pressure 
% (based on power production)

for iVar = 1:nVar
    iVar1 = find(NumWECs == NumWECs_mesh(iVar));
    iVar2 = find(VperWEC == VperWEC_mesh(iVar));
    iVar3 = find(p_nom == p_nom_mesh(iVar));

    PP_WEC(iVar1,iVar2,iVar3,SS) = PP_WEC_array(iVar,SS);
    PP_wp(iVar1,iVar2,iVar3,SS) = PP_wp_array(iVar,SS);
    PP_mLoss(iVar1,iVar2,iVar3,SS) = PP_mLoss_array(iVar,SS);
    PP_gen(iVar1,iVar2,iVar3,SS) = PP_gen_array(iVar,SS);
    PP_hPRV(iVar1,iVar2,iVar3,SS) = PP_hPRV_array(iVar,SS);
    dpdt_max(iVar1,iVar2,iVar3,SS) = dpdt_max_array(iVar,SS);

end

for iVar1 = 1:nVar1
    for iVar2 = 1:nVar2
        iVar3 = ...
         find(PP_gen(iVar1,iVar2,:,SS) == max(PP_gen(iVar1,iVar2,:,SS)),1);
        if isnan(max(PP_gen(iVar1,iVar2,:,SS)))
            iVar3 = 1;
        end

        PP_WEC_opt(iVar1,iVar2,SS) = PP_WEC(iVar1,iVar2,iVar3,SS);
        PP_wp_opt(iVar1,iVar2,SS) = PP_wp(iVar1,iVar2,iVar3,SS);
        PP_mLoss_opt(iVar1,iVar2,SS) = PP_mLoss(iVar1,iVar2,iVar3,SS);
        PP_gen_opt(iVar1,iVar2,SS) = PP_gen(iVar1,iVar2,iVar3,SS);
        PP_hPRV_opt(iVar1,iVar2,SS) = PP_hPRV(iVar1,iVar2,iVar3,SS);
        dpdt_max_opt(iVar1,iVar2,SS) = dpdt_max(iVar1,iVar2,iVar3,SS);

    end
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




NumWECs_to_plot = [1 2 3 4 5];
I = [];
for i = 1:numel(NumWECs_to_plot)
    I(i) = find(NumWECs == NumWECs_to_plot(i));
end

p = semilogx([0 0], -[99 99]);
p.HandleVisibility = 'off';

hold on
for i = 1:numel(I)
    scatter(VperWEC,PP_gen_opt(I(i),:,SS)./NumWECs_to_plot(i)/PP_max,100,color(i,:),'Marker','x','LineWidth',2)
    hold on
    legStr(i) = {[num2str(NumWECs(I(i))),' WECs']};
end
yLim = ylim;
ylim([0,yLim(2)])
ylim([0,1])

ylabel('Elec. power, mean (kW)', ...
    'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)
xlabel('Volume (1000L)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)


grid on

PowerFactor = PP_max/(D_m_base*par.motor.w_max*(maxPressure-par.control.p_l_nom))

title(['Mean Power Production Vs. ', ...
    'Accumulator Volume per WEC:',newline, ...
    'Power Factor of ',num2str(PowerFactor,3)],...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname',fontName)

ax = gca;
ax.FontName = fontName;
ax.FontSize = axFontSize;

leg = legend(legStr);
leg.FontSize = axFontSize;
leg.FontName = fontName;
set(leg, 'Location', 'best')