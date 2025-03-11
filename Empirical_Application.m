clear; clc; rng(100)
addpath('functions','data')

print_figures = 1; % set to 1 if figures are to be printed into the output folder

%% Load Oil Market Dataset
opts = detectImportOptions("dataset_oilmarket.xlsx");
opts.VariableTypes(2:end)={'double'};
opts.Sheet = 'Export';
opts.VariableNames(1)={'dates'};
data = readtable("dataset_oilmarket.xlsx",opts);

%% Load G17 dataset: industrial production indices and relative importance weights
opts = detectImportOptions("dataset_IP.csv");
opts.VariableTypes(2:end)={'double'};
opts.VariableTypes(1) = {'char'};
opts.VariableNames(1)={'date'};
data_IP = readtable("dataset_IP.csv",opts);
data_IP.date = datetime(data_IP.date)+calmonths(1)-caldays(1);
data_IP = table2timetable(data_IP);
data_IP = data_IP(data_IP.date<=datetime(2023,12,31),:);
is_IP = contains(data_IP.Properties.VariableNames,'IP_');
N = sum(is_IP);
is_RIW = contains(data_IP.Properties.VariableNames,'RIW_');
varnames_IP  = data_IP.Properties.VariableNames(is_IP);
varnames_RIW = data_IP.Properties.VariableNames(is_RIW);
naics_codes_IP = NaN(N,1);
for i = 1:size(varnames_IP,2)
    varname_plus_desc = split(string(varnames_IP{i}),"_S");
    if ~isempty(char(varname_plus_desc(2)))
        varnames_IP_desc{i}  = char(varname_plus_desc(2));
    else
        varnames_IP_desc{i}  = char(varname_plus_desc(3));
    end
    varnames_IP{i} = char(strcat(varname_plus_desc(1),"_S"));
    naics_codes_IP(i) = str2double( regexp(varname_plus_desc(1),'\d*','Match'));
end
idx_MFG = and(naics_codes_IP>=311,naics_codes_IP<=339);
idx_Mining =  naics_codes_IP>=2111;
N_mfg = sum(idx_MFG);  N_mining = sum(idx_Mining);
varnames_IP_mfg = varnames_IP(idx_MFG); varnames_IP_mining = varnames_IP(idx_Mining);
varnames_RIW_mfg = varnames_RIW(idx_MFG); varnames_RIW_mining = varnames_RIW(idx_Mining);
Xall_IP = table2array(data_IP(:,is_IP));
Rall = table2array(data_IP(:,is_RIW));
Xall_mfg = Xall_IP(:,idx_MFG); Xall_mining = Xall_IP(:,idx_Mining);
Rall_mfg = Rall(:,idx_MFG); Rall_mining = Rall(:,idx_Mining);
Rall_mfg = Rall_mfg./sum(Rall_mfg,2); Rall_mining = Rall_mining./sum(Rall_mining,2);
idx_est_indstry =  and(data_IP.date>=datetime(1974,1,1),data_IP.date<= datetime(2023,12,31) );

%% Colors for Figures 
color_1 = [232,232,232]./255;
color_2 = [128,128,128]./255; 
facealpha = 0.4;
facealpha2 = 0.25;


%% Find a Bandwidth on Pre-Covid dataset
date_bs = data.dates;
IV_DK = data.DK_IV;
IV_DK(isnan(IV_DK))=0;
X_all = [IV_DK, 100*log(data.WTI./data.USCPI), 100*log(data.WorldOilProduction_mbpd_),  ...
    100*log(data.CrudeStocksProxySA),  100*log(data.WorldIP) , 100*log(data.USMFGIP), 100*log(data.USMiningIP) ];
varnames = {'RPO','Crude production', 'Crude Stocks','World IP', 'US Manufacturing IP', 'US Mining IP'};%,'US IP (MFG)','US IP (Petroleum and Coal')
idx_est_CV = and(date_bs>=datetime(1974,1,1),date_bs<= datetime(2019,12,31) ) ; % Pre-Covid dataset
X_all_CV = X_all(idx_est_CV,:);
Hgrid = length(X_all_CV).^(0.5:0.005:0.9);
CV_var_choice = 2:size(X_all_CV,2) ; % Which variables to include in the conditional forecast exercise
CV_fraction_sample = 1/2; % Where to start in the sample
nthin = 1; hfc = 1; % One Step Ahead Forecast
[Hest, CVs] = choose_H_CF( [X_all_CV(:,1),X_all_CV(:,2:end)] , 13, 1, Hgrid, nthin, hfc, CV_var_choice, CV_fraction_sample);



%% Replicate Figure 7 (baseline results), but choosing a slightly lower bandwith of H  = 150 (see text)
H = 150;
idx_est =  and(date_bs>=datetime(1974,1,1),date_bs<= datetime(2023,12,31) );
date_bs = date_bs(idx_est);
X_all = X_all(idx_est,:);
p = 13; h = 5*12;
ngrid = 6;
alpha = 0.1; alpha2 = 0.32;
% Create COVID dummies:
date_U = date_bs(p+1:end);
T = length(date_bs);
idx =  ceil(linspace(ceil(T/20), ceil( (T-p-1)/20*19 ), ngrid) );
dates_est = date_U(idx);
dates_dummies = and(date_bs>datetime(2020,2,1),date_bs<datetime(2022,1,1)+calmonths(12)); % plus calmonths(p)???+calmonths(VAR.p)
X_exo = zeros(T,1+sum(dates_dummies));
X_exo(:,1)=1;
X_exo(dates_dummies,2:end) = eye(sum(dates_dummies));
% Estimate model (for two alphas)
[ output_externalIV_IP ]   = tvp_svar_iv_kernel(X_all(:,2:end),  X_all(:,1) , p, H , alpha , X_exo, idx , h );
[ output_externalIV_IP_2 ] = tvp_svar_iv_kernel(X_all(:,2:end),  X_all(:,1) , p, H , alpha2 , X_exo, idx , h );

% Print the Portmanteau Test
disp('Portmanteau Test: dates')
disp(dates_est')
disp('Portmanteau Test: Statistic')
disp(round(output_externalIV_IP.Portmanteau_test_testat',2) )
disp('Portmanteau Test: DoF')
disp(output_externalIV_IP.Portmanteau_test_dof')
disp('Portmanteau Test: p-val')
disp(output_externalIV_IP.Portmanteau_test_pval' )


% Print the F-Tests of invertibility
disp('Wald test for invertibility: dates')
disp(dates_est')
disp('Wald test for invertibility: test statistics')
disp(round(output_externalIV_IP.Waldtest',2))
disp('Wald test for invertibility: p-values')
disp(round(output_externalIV_IP.Waldtest_pval',2))
disp('F test for invertibility: test statistics')
disp(round(output_externalIV_IP.Ftest',2))
disp('F test for invertibility: p-values')
disp(round(output_externalIV_IP.Ftest_pval_1',2))
disp(round(output_externalIV_IP.Ftest_pval_2',2))



% Plot of IRFs (main estimates)
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
        p1 = plot((0:h)', output_externalIV_IP.IVSVAR(:,ii,i), 'k','Linewidth',2,'LineStyle','-');
        p2 = plot((0:h)', [output_externalIV_IP.IVSVAR_lb_osw(:,ii,i),output_externalIV_IP.IVSVAR_ub_osw(:,ii,i)], 'k--','Linewidth',1);
        Y_2 = [output_externalIV_IP.IVSVAR_lb_osw(:,ii,i)', fliplr( output_externalIV_IP.IVSVAR_ub_osw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );
        Y_3 = [output_externalIV_IP_2.IVSVAR_lb_osw(:,ii,i)', fliplr( output_externalIV_IP_2.IVSVAR_ub_osw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','--');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );
        ylim([min(min(output_externalIV_IP.IVSVAR_lb_osw(:,ii,:))),max(max(output_externalIV_IP.IVSVAR_ub_osw(:,ii,:)))])
        title(varnames(ii),'interpreter','latex')
        p3 =  plot((0:h)', output_externalIV_IP.IVSVAR_cp(:,ii), 'k','Linewidth',2,'LineStyle',':');
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
if print_figures == 1
    print(hfig1, 'output/IRFs_baseline', '-dpdf')
end
 
% Plot of alpha_t
dates_plot = date_U(idx);
dates_plot.Format = 'MM/yy';
xgraph_alpha = [dates_plot', fliplr(dates_plot')];

hfig2 = figure(2);
grid on; hold on;
plot(dates_plot,output_externalIV_IP.alpha,'k','linewidth',2)
ylim([0,.5])
Y_2 =  [output_externalIV_IP.alpha_lb', fliplr(output_externalIV_IP.alpha_ub')];
bounds = fill(xgraph_alpha,Y_2,color_2,'LineStyle','none');
set(bounds,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );
Y_3 = [output_externalIV_IP_2.alpha_lb', fliplr( output_externalIV_IP_2.alpha_ub')  ];
bounds2 = fill(xgraph_alpha,Y_3,color_2,'LineStyle','--');
set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );

title('Variation in $\alpha_t$','Interpreter','latex')
set(hfig2,'PaperPositionMode','auto')
set(hfig2, 'Position', [30 50 600 200])
hfig2 = tightfig(hfig2);
if print_figures == 1
    print(hfig2,'output\figure_alpha', '-dpdf')
end


%% Robustness exercises
% Internal instrument VAR (relative IRFs)
idx_shock = 1;
[ output_internalIV_IP ] = tvp_svar_internal_kernel(X_all(:,2:end),  X_all(:,1) , p, H , alpha , X_exo, idx , h, idx_shock );
[ output_internalIV_IP_2 ] = tvp_svar_internal_kernel(X_all(:,2:end),  X_all(:,1) , p, H , alpha2 , X_exo, idx , h, idx_shock );
shocksize_internal_IV = output_externalIV_IP.IVSVAR(1,1,output_internalIV_IP.t_stand_idx); % This scalar is used to rescale all the IRFs.
% External IV: different bandwidths bandwidths
H2 = 110; % Roughly the estimated Bandwith
[ output_externalIV_IP_H2 ] = tvp_svar_iv_kernel(X_all(:,2:end),  X_all(:,1) , p, H2 , alpha , X_exo, idx , h );
H3 = 190;
[ output_externalIV_IP_H3 ] = tvp_svar_iv_kernel(X_all(:,2:end),  X_all(:,1) , p, H3 , alpha , X_exo, idx , h );


% Plot
hfig3 = figure(3);
a = 1;
for ii = size(varnames,2)-1:size(varnames,2)
    for i = 1:size(dates_est,1)
        subplot(4,size(dates_est,1),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', shocksize_internal_IV.*output_internalIV_IP.VARpm(:,ii,i), 'k','Linewidth',2);
        %p2 = plot((0:h)', shocksize_internal_IV.*[output_internalIV_IP.VARpm_lb_msw(:,ii,i),output_internalIV_IP.VARpm_ub_msw(:,ii,i)], 'k--','Linewidth',1);
        Y_2 = shocksize_internal_IV.*[output_internalIV_IP.VARpm_lb_msw(:,ii,i)', fliplr( output_internalIV_IP.VARpm_ub_msw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );
        Y_3 = shocksize_internal_IV.*[output_internalIV_IP_2.VARpm_lb_msw(:,ii,i)', fliplr( output_internalIV_IP_2.VARpm_ub_msw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','--');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );
        ylim([min(min(shocksize_internal_IV.*output_internalIV_IP.VARpm_lb_msw(:,ii,:))),max(max(shocksize_internal_IV.*output_internalIV_IP.VARpm_ub_msw(:,ii,:)))])
        p3 = plot((0:h)', output_externalIV_IP.IVSVAR(:,ii,i), 'k--','Linewidth',2);
        title(varnames(ii),'interpreter','latex')
        ylabel(datestr(dates_est(i),'yyyy'))
        if i == 1
            legend([p1,p3],{'Internal IV VAR','IV-SVAR'},'interpreter','latex','location','northwest')
        end
        a = a + 1;
    end
end
for ii = size(varnames,2)-1:size(varnames,2)
    for i = 1:size(dates_est,1)
        subplot(4,size(dates_est,1),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', output_externalIV_IP.IVSVAR(:,ii,i), 'k','Linewidth',2);
       % p2 = plot((0:h)', [output_externalIV_IP.IVSVAR_lb_osw(:,ii,i),output_externalIV_IP.IVSVAR_ub_osw(:,ii,i)], 'k--','Linewidth',1);
        Y_2 = [output_externalIV_IP.IVSVAR_lb_osw(:,ii,i)', fliplr( output_externalIV_IP.IVSVAR_ub_osw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );
        Y_3 = [output_externalIV_IP_2.IVSVAR_lb_osw(:,ii,i)', fliplr( output_externalIV_IP_2.IVSVAR_ub_osw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','--');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );
        ylim([min(min(output_externalIV_IP.IVSVAR_lb_osw(:,ii,:))),max(max(output_externalIV_IP.IVSVAR_ub_osw(:,ii,:)))])
        title(varnames(ii),'interpreter','latex')
        ylabel(datestr(dates_est(i),'yyyy'))
        p3 = plot((0:h)', output_externalIV_IP_H2.IVSVAR(:,ii,i), 'k','LineStyle',':','Linewidth',2);
        p4 = plot((0:h)', output_externalIV_IP_H3.IVSVAR(:,ii,i), 'k','LineStyle','--','Linewidth',2);
        if i == 1
            legend([p1,p3,p4],{'IV-SVAR ($H=150$)','IV-SVAR ($H=110$)','IV-SVAR ($H=190$)'},'interpreter','latex','location','northwest')
        end
        a = a + 1;
    end
end
set(gcf,'PaperPositionMode','auto')
set(hfig3, 'Position', [30 50 1300 900])
hfig3 = tightfig(hfig3) ;
if print_figures == 1
    print(hfig3, 'output/IRFs_baseline_robustness', '-dpdf')
end


% Print all the IRFs of the internal instrument VAR
xgraph = [(0:h), fliplr((0:h))];
hfig100 = figure(100);
a = 1;
for i = 1:size(dates_est,1)
    for ii = 1:size(varnames,2)
        subplot(size(dates_est,1),size(varnames,2),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', shocksize_internal_IV.*output_internalIV_IP.VARpm(:,ii,i), 'k','Linewidth',2);
        p2 = plot((0:h)', shocksize_internal_IV.*[output_internalIV_IP.VARpm_lb_msw(:,ii,i),output_internalIV_IP.VARpm_ub_msw(:,ii,i)], 'k--','Linewidth',1);
        Y_2 = [shocksize_internal_IV.*output_internalIV_IP.VARpm_lb_msw(:,ii,i)', fliplr( shocksize_internal_IV.*output_internalIV_IP.VARpm_ub_msw(:,ii,i)')  ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );
        ylim([min(min(shocksize_internal_IV.*output_internalIV_IP.VARpm_lb_msw(:,ii,:))),max(max(shocksize_internal_IV.*output_internalIV_IP.VARpm_ub_msw(:,ii,:)))])
        title(varnames(ii),'interpreter','latex')
        p3 =  plot((0:h)', shocksize_internal_IV.*output_internalIV_IP.VARpm_cp(:,ii), 'k:','Linewidth',2);
       % p3b =  plot((0:h)', shocksize_internal_IV.*[output_internalIV_IP.VARpm_cp_ub_msw(:,ii),output_internalIV_IP.VARpm_cp_lb_msw(:,ii)], 'r--','Linewidth',1);
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
set(hfig100, 'Position', [30 50 800 800])
hfig100 = tightfig(hfig100) ;
if print_figures == 1
    print(hfig100, 'output/IRFs_baseline_internal_IV', '-dpdf')
end




%% Industry-level exercises  
outputs_all = {};
outputs_all_2 = {};
for i = 1:N
    atom_prefix = char(varnames_IP(i));
    atom = atom_prefix(4:end);
    if idx_MFG(i) == 1
        id_IP  = find(contains(varnames_IP_mfg,atom));
        id_RIW = find(contains(varnames_RIW_mfg,atom));
        IP_i = Xall_mfg(:,id_IP);
        % Manufacturing industry
        IP_exi = 100*compute_special_IP(Xall_mfg(:,setdiff(1:N_mfg,id_IP)),Rall_mfg(:,setdiff(1:N_mfg,id_RIW)));
        X_IP = [X_all(:,end),100*log(IP_exi(idx_est_indstry)), 100*log(IP_i(idx_est_indstry)) ];
    elseif idx_MFG(i)==0
        id_IP  = find(contains(varnames_IP_mining,atom));
        id_RIW = find(contains(varnames_RIW_mining,atom));
        IP_i = Xall_mining(:,id_IP);
        % Mining industry
        IP_exi = 100*compute_special_IP(Xall_mining(:,setdiff(1:N_mining,id_IP)),Rall_mining(:,setdiff(1:N_mining,id_RIW)));
        X_IP = [X_all(:,end-1), 100*log(IP_exi(idx_est_indstry)), 100*log(IP_i(idx_est_indstry)) ];
    end

    % Match with macro data
    X_all_i = [ X_all(:,1:end-2), X_IP ];
    [ output_externalIV ] = tvp_svar_iv_kernel_fast(X_all_i(:,2:end),  X_all_i(:,1) , p, H , alpha , X_exo, idx , h );
    [ output_externalIV_2 ] = tvp_svar_iv_kernel_fast(X_all_i(:,2:end),  X_all_i(:,1) , p, H , alpha2 , X_exo, idx , h );

    %  Save IRFs for each industry:
    outputs_all{i} = output_externalIV;
    outputs_all_2{i} = output_externalIV_2;
end


%% Make figures by aggregates
naics_codes = readtable('data/2-6 digit_2022_Codes.csv'); 
MainPaper = [2111, 324,  325, 326 , 336, 331, 332 ];
a = 1;
hfig_main_paper = figure(4); 
for i =  [3,6]
    for ii = 1:size(MainPaper,2)
        subplot(2,size(MainPaper,2),a)
        code_i = MainPaper(ii);
        id_i = find(naics_codes_IP==code_i);
        output_i = outputs_all{id_i};
        hold on; grid on;
        p1 = plot((0:h)', output_i.IVSVAR(:,end,i), 'k','Linewidth',2);
        Y_2 = [output_i.IVSVAR_lb_osw(:,end,i)', fliplr( output_i.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );
        p2 = plot((0:h)', [output_i.IVSVAR_lb_osw(:,end,i),output_i.IVSVAR_ub_osw(:,end,i)], 'k--','Linewidth',1);
        output_i_2 = outputs_all_2{id_i};
        Y_3 = [output_i_2.IVSVAR_lb_osw(:,end,i)', fliplr( output_i_2.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );
        title_i =  erase(naics_codes.x2022NAICSUSTitle(find(naics_codes.x2022NAICSUSCode==naics_codes_IP(id_i))),'Manufacturing');
        title(title_i,'interpreter','latex')
        p3 =  plot((0:h)', output_i.IVSVAR_cp(:,end), 'k:','Linewidth',2);
        xlim([0,h]) 
        if mod(a,size(MainPaper,2))==1
            ylabel(datestr(dates_est(i),'yyyy'))
        end
        if a == 1
            legend([p1,p3],{'TVP','CP'},'interpreter','latex')
        end
        a = a + 1;
    end
end
set(gcf,'PaperPositionMode','auto')
set(hfig_main_paper, 'Position', [30 50 1400 400])
hfig_main_paper = tightfig(hfig_main_paper) ;
if print_figures == 1
    print(hfig_main_paper, 'output/IRFs_by_industry_maintext', '-dpdf')  
    %print(hfig_main_paper, 'output/IRFs_by_industry_maintext',  '-dpng','-r500')
end

% Durables
Durables1 = [321, 327, 331, 332 , 333  ];
a = 1;
hfig_d1  = figure(5);
for i = 1:size(dates_est,1)
    for ii = 1:size(Durables1,2)
        subplot(size(dates_est,1),size(Durables1,2),a)
        code_i = Durables1(ii);
        id_i = find(naics_codes_IP==code_i);
        output_i = outputs_all{id_i};
        hold on; grid on;
        p1 = plot((0:h)', output_i.IVSVAR(:,end,i), 'k','Linewidth',2);
        Y_2 = [output_i.IVSVAR_lb_osw(:,end,i)', fliplr( output_i.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );

        p2 = plot((0:h)', [output_i.IVSVAR_lb_osw(:,end,i),output_i.IVSVAR_ub_osw(:,end,i)], 'k--','Linewidth',1);
        output_i_2 = outputs_all_2{id_i};
        Y_3 = [output_i_2.IVSVAR_lb_osw(:,end,i)', fliplr( output_i_2.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );
        title_i =  erase(naics_codes.x2022NAICSUSTitle(find(naics_codes.x2022NAICSUSCode==naics_codes_IP(id_i))),'Manufacturing');
        title(title_i,'interpreter','latex')
        p3 =  plot((0:h)', output_i.IVSVAR_cp(:,end), 'k:','Linewidth',2);
        xlim([0,h])
        if mod(a,size(Durables1,2))==1
            ylabel(datestr(dates_est(i),'yyyy'))
        end
        if a == 1
            legend([p1,p3],{'TVP','CP'},'interpreter','latex')
        end
        a = a + 1;
    end
end
set(gcf,'PaperPositionMode','auto')
set(hfig_d1, 'Position', [30 50 1400 800])
hfig_d1 = tightfig(hfig_d1) ;
if print_figures == 1
    print(hfig_d1, 'output/IRFs_durables1', '-dpdf')
    %print(hfig_d1, 'output/IRFs_durables1', '-dpng','-r500')  
end

Durables2 = [334, 335, 336, 337, 339 ];
a = 1;
hfig_d2  = figure(6);
for i = 1:size(dates_est,1)
    for ii = 1:size(Durables2,2)
        subplot(size(dates_est,1),size(Durables2,2),a)
        code_i = Durables2(ii);
        id_i = find(naics_codes_IP==code_i);
        output_i = outputs_all{id_i};
        hold on; grid on;
        p1 = plot((0:h)', output_i.IVSVAR(:,end,i), 'k','Linewidth',2);
        Y_2 = [output_i.IVSVAR_lb_osw(:,end,i)', fliplr( output_i.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );

        p2 = plot((0:h)', [output_i.IVSVAR_lb_osw(:,end,i),output_i.IVSVAR_ub_osw(:,end,i)], 'k--','Linewidth',1);
        output_i_2 = outputs_all_2{id_i};
        Y_3 = [output_i_2.IVSVAR_lb_osw(:,end,i)', fliplr( output_i_2.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );
        title_i =  erase(naics_codes.x2022NAICSUSTitle(find(naics_codes.x2022NAICSUSCode==naics_codes_IP(id_i))),'Manufacturing');
        title(title_i,'interpreter','latex') 
        p3 =  plot((0:h)', output_i.IVSVAR_cp(:,end), 'k:','Linewidth',2);
        xlim([0,h]) 
        if mod(a,size(Durables2,2))==1
            ylabel(datestr(dates_est(i),'yyyy'))
        end
        if a == 1
            legend([p1,p3],{'TVP','CP'},'interpreter','latex')
        end
        a = a + 1;
    end
end
set(gcf,'PaperPositionMode','auto')
set(hfig_d2, 'Position', [30 50 1400 800])
hfig_d2 = tightfig(hfig_d2) ;
if print_figures == 1
    print(hfig_d2, 'output/IRFs_durables2', '-dpdf')
    %print(hfig_d2, 'output/IRFs_durables2', '-dpng','-r500')   
end


% Nondurables
NonDurables1 = [311, 312, 313, 314, 315 ];
a = 1;
hfig_nd = figure(7);
for i = 1:size(dates_est,1)
    for ii = 1:size(NonDurables1,2)
        subplot(size(dates_est,1),size(NonDurables1,2),a)
        code_i = NonDurables1(ii);
        id_i = find(naics_codes_IP==code_i);
        output_i = outputs_all{id_i};
        hold on; grid on;
        p1 = plot((0:h)', output_i.IVSVAR(:,end,i), 'k','Linewidth',2);
        Y_2 = [output_i.IVSVAR_lb_osw(:,end,i)', fliplr( output_i.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );

        p2 = plot((0:h)', [output_i.IVSVAR_lb_osw(:,end,i),output_i.IVSVAR_ub_osw(:,end,i)], 'k--','Linewidth',1);
        output_i_2 = outputs_all_2{id_i};
        Y_3 = [output_i_2.IVSVAR_lb_osw(:,end,i)', fliplr( output_i_2.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );

        title_i =  erase(naics_codes.x2022NAICSUSTitle(find(naics_codes.x2022NAICSUSCode==naics_codes_IP(id_i))),'Manufacturing');
        title(title_i,'interpreter','latex')
        p3 =  plot((0:h)', output_i.IVSVAR_cp(:,end), 'k:','Linewidth',2);
        xlim([0,h])
        if mod(a,size(NonDurables1,2))==1
            ylabel(datestr(dates_est(i),'yyyy'))
        end

        if a == 1
            legend([p1,p3],{'TVP','CP'},'interpreter','latex')
        end
            
        a = a + 1;
    end
end

set(gcf,'PaperPositionMode','auto')
set(hfig_nd, 'Position', [30 50 1400 800])
hfig_nd = tightfig(hfig_nd) ;
if print_figures == 1
    print(hfig_nd, 'output/IRFs_nondurables1', '-dpdf')
    %print(hfig_nd, 'output/IRFs_nondurables1', '-dpng','-r500')  
end

NonDurables2 = [316, 322,323, 324, 325, 326];
a = 1;
hfig_nd2 = figure(8);
for i = 1:size(dates_est,1)
    for ii = 1:size(NonDurables2,2)
        subplot(size(dates_est,1),size(NonDurables2,2),a)
        code_i = NonDurables2(ii);
        id_i = find(naics_codes_IP==code_i);
        output_i = outputs_all{id_i};
        hold on; grid on;
        p1 = plot((0:h)', output_i.IVSVAR(:,end,i), 'k','Linewidth',2);
        Y_2 = [output_i.IVSVAR_lb_osw(:,end,i)', fliplr( output_i.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );

        p2 = plot((0:h)', [output_i.IVSVAR_lb_osw(:,end,i),output_i.IVSVAR_ub_osw(:,end,i)], 'k--','Linewidth',1);
        output_i_2 = outputs_all_2{id_i};
        Y_3 = [output_i_2.IVSVAR_lb_osw(:,end,i)', fliplr( output_i_2.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );

        title_i =  erase(naics_codes.x2022NAICSUSTitle(find(naics_codes.x2022NAICSUSCode==naics_codes_IP(id_i))),'Manufacturing');
        title(title_i,'interpreter','latex')
        p3 =  plot((0:h)', output_i.IVSVAR_cp(:,end), 'k:','Linewidth',2);
        xlim([0,h])
        if mod(a,size(NonDurables2,2))==1
            ylabel(datestr(dates_est(i),'yyyy'))
        end

        if a == 1
            legend([p1,p3],{'TVP','CP'},'interpreter','latex')
        end
        a = a + 1;
    end
end
set(gcf,'PaperPositionMode','auto')
set(hfig_nd2, 'Position', [30 50 1400 800])
hfig_nd2 = tightfig(hfig_nd2) ;
if print_figures == 1
    print(hfig_nd2, 'output/IRFs_nondurables2', '-dpdf')
    %print(hfig_nd2, 'output/IRFs_nondurables2', '-dpng','-r500')
end

% Energy

Energy = [2111, 2122,2123,2131, 2121];
a = 1;
hfig_mining = figure(9);
for i = 1:size(dates_est,1)
    for ii = 1:size(Energy,2)
        subplot(size(dates_est,1),size(Energy,2),a)
        code_i = Energy(ii);
        id_i = find(naics_codes_IP==code_i);
        output_i = outputs_all{id_i};
        hold on; grid on;
        p1 = plot((0:h)', output_i.IVSVAR(:,end,i), 'k','Linewidth',2);
        Y_2 = [output_i.IVSVAR_lb_osw(:,end,i)', fliplr( output_i.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 );

        p2 = plot((0:h)', [output_i.IVSVAR_lb_osw(:,end,i),output_i.IVSVAR_ub_osw(:,end,i)], 'k--','Linewidth',1);
        output_i_2 = outputs_all_2{id_i};
        Y_3 = [output_i_2.IVSVAR_lb_osw(:,end,i)', fliplr( output_i_2.IVSVAR_ub_osw(:,end,i)') ];
        bounds2 = fill(xgraph,Y_3,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha ,'EdgeAlpha',facealpha );

        title_i =  erase(naics_codes.x2022NAICSUSTitle(find(naics_codes.x2022NAICSUSCode==naics_codes_IP(id_i))),'Manufacturing');
        title(title_i,'interpreter','latex')
        p3 =  plot((0:h)', output_i.IVSVAR_cp(:,end), 'k:','Linewidth',2);
        xlim([0,h])
        if mod(a,size(Energy,2))==1
            ylabel(datestr(dates_est(i),'yyyy'))
        end

        if a == 1
            legend([p1,p3],{'TVP','CP'},'interpreter','latex')
        end

        a = a + 1;
    end
end
set(gcf,'PaperPositionMode','auto')
set(hfig_mining, 'Position', [30 50 1400 800])
hfig_mining = tightfig(hfig_mining) ;
if print_figures == 1
    print(hfig_mining, 'output/IRFs_mining', '-dpdf')
    %print(hfig_mining, 'output/IRFs_mining', '-dpng','-r500')
end











