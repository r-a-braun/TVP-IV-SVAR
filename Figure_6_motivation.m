clear; clc;
addpath('functions','data')

[tempData, tempText] = xlsread('Table31oil.xlsx','Annual Data'); 
T = array2table(tempData,"VariableNames",tempText(11,1:end));
T.("Annual Average")=datetime( T.("Annual Average"),12,31);
T = table2timetable(T); 
%% Figure 
hfig = figure(1);
hold on; grid on;
production = T.("Total Petroleum Field Production") + T.("Biofuels Plant Net Production") + T.("Petroleum Processing Gain");
p1 = plot(T.("Annual Average"),T.("Petroleum Products Supplied")/1000,'k','LineWidth',2,'LineStyle','-');
p2 = plot(T.("Annual Average"),T.("Petroleum Production")/1000,'k','LineWidth',2,'LineStyle','--');
p3 = plot(T.("Annual Average"),T.("Petroleum Net Imports")/1000,'k','LineWidth',2,'LineStyle',':');
ylabel('Million barrels per day') 
%recessionplot 
legend([p1,p2,p3],{'Consumption','Production','Net Imports'},'interpreter','latex','Location','northwest') 
set(hfig, 'Position', [30 50 600 400])
hfig = tightfig(hfig) ;
print(hfig,  'output/us_petroleum_developments' , '-dpdf')
close(hfig)




