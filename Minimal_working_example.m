% This code illustrates the use of the TVP IV-SVAR estimator 
% by replicating only Figure 7
clear; clc; rng(100)
addpath('functions','data') 

%% Load Oil Market Dataset
opts = detectImportOptions("dataset_oilmarket.xlsx");
opts.VariableTypes(2:end)={'double'};
opts.Sheet = 'Export';
opts.VariableNames(1)={'dates'};
data = readtable("dataset_oilmarket.xlsx",opts);
date_sample = data.dates;
IV_DK = data.DK_IV;
IV_DK(isnan(IV_DK))=0;
X_all = [IV_DK, 100*log(data.WTI./data.USCPI), 100*log(data.WorldOilProduction_mbpd_),  ...
    100*log(data.CrudeStocksProxySA),  100*log(data.WorldIP) , 100*log(data.USMFGIP), 100*log(data.USMiningIP) ];
varnames = {'RPO','Crude production', 'Crude Stocks','World IP', 'US Manufacturing IP', 'US Mining IP'};%,'US IP (MFG)','US IP (Petroleum and Coal')

%% Choose a bandwidth (H), a sample (idx_est), lag length (p), IRF horizons (h), confidence level (alpha), 
% a grid of periods where the model is estimated
% and set up COVID dummies.
H = 150;
idx_est =  and(date_sample>=datetime(1974,1,1),date_sample<= datetime(2023,12,31) );
date_sample = date_sample(idx_est);
X_all = X_all(idx_est,:);
p = 13; h = 5*12;
ngrid = 6;
alpha = 0.1;  
T = length(date_sample);
idx =  ceil(linspace(ceil(T/20), ceil( (T-p-1)/20*19 ), ngrid) ); 
date_U = date_sample(p+1:end);
dates_est = date_U(idx);
dates_dummies = and(date_sample>datetime(2020,2,1),date_sample<datetime(2022,1,1)+calmonths(12));  
X_exo = zeros(T,1+sum(dates_dummies));
X_exo(:,1) = 1;
X_exo(dates_dummies,2:end) = eye(sum(dates_dummies));

%% Estimate the model 
[ output_IVSVAR ]   = tvp_svar_iv_kernel(X_all(:,2:end),  X_all(:,1) , p, H , alpha , X_exo, idx , h ); 
 
%% Plot the output 
colors = linspecer(3+4); 
color_2 = colors(2,:); 
facealpha = 0.4;
xgraph = [(0:h), fliplr((0:h))];
hfig1 = figure(1);
a = 1;
T_figure_all = table((0:h)');
T_figure_all.Properties.VariableNames = {'horizon'};
for i = 1:size(dates_est,1)
    for ii = 1:size(varnames,2)
        subplot(size(dates_est,1),size(varnames,2),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', output_IVSVAR.IVSVAR(:,ii,i), 'color',color_2,'Linewidth',2);
        p2 = plot((0:h)', [output_IVSVAR.IVSVAR_lb_osw(:,ii,i),output_IVSVAR.IVSVAR_ub_osw(:,ii,i)], 'k--','Linewidth',1);
        Y_2 = [output_IVSVAR.IVSVAR_lb_osw(:,ii,i)', fliplr( output_IVSVAR.IVSVAR_ub_osw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha ); 
        ylim([min(min(output_IVSVAR.IVSVAR_lb_osw(:,ii,:))),max(max(output_IVSVAR.IVSVAR_ub_osw(:,ii,:)))])
        title(varnames(ii),'interpreter','latex')
        p3 =  plot((0:h)', output_IVSVAR.IVSVAR_cp(:,ii), 'r','Linewidth',2);
        if a == 1
            legend([p1,p3],{'TVP','CP'},'interpreter','latex')
        end
        xlim([0,h])
        if mod(a,size(varnames,2))==1
            ylabel(datestr(dates_est(i),'yyyy'))
        end
        a = a + 1;
    end
end
set(gcf,'PaperPositionMode','auto')
set(hfig1, 'Position', [30 50 800 800])
hfig1 = tightfig(hfig1) ;
