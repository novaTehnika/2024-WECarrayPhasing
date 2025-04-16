black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

linestyles = {'-', '--', ':', '-.'};

fontName = 'Palatino Linotype';
supTitleFontSize = 10;
subTitleFontSize = 10;
axFontSize = 8;
legFontSize = 8;
bottomEdge = 1*2.5;
leftEdge = 3*2.5;
width = 6; % one column: 3+9/16, two column: 7.5
height = 5;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'centimeters';
fig.Position = [leftEdge bottomEdge width height ];
set(fig,'defaultAxesColorOrder',[black; black]);

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = fontName;
ax1.FontSize = axFontSize;

SS = 2;
p = plot(T_c_data(SS,:)*1e-6,PP_w_data(SS,:)*1e-3,'Color',black,'LineWidth',lineWidth);


ylabel('Power (kW)', ...
    'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)
xlabel('Coulomb damping torque (MNm)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname',fontName)


grid on

% title(['Mean power production per WEC', ...
%     '',newline, ...
%     num2str(NumWECs), ' WECs',...
%     ' and ','a power factor of ',num2str(PowerFactor,3)],...
% 'Interpreter','latex','FontSize',supTitleFontSize,'fontname',fontName)

% title([num2str(NumWECs), ' WECs',...
%     ' and ','a power factor of ',num2str(PowerFactor,3)],...
% 'Interpreter','latex','FontSize',supTitleFontSize,'fontname',fontName)

ax = gca;
ax.FontName = fontName;
ax.FontSize = axFontSize;

% leg = legend(legStr);
% leg.FontSize = legFontSize;
% leg.FontName = fontName;
% set(leg, 'Location', 'best')