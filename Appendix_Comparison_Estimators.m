clear; clc;
rng(06112023)
%% This file runs the exercise to compare the kernel estimator to a path estimator (Mueller and Petalas type) and the Bayesian estimator 3
%% See Appendix F
T = 500; n = 2; p =1;
A1 = [.8,-.05;.2,.7];
A2 = A1; 
B1 = [1,.2;1,1];
B2 = [.5,.4;-1,1]; 
phi1 = 0.2; phi2 = 0.1;
sig_v = 0.1;
y = zeros(T,n);
z = zeros(T,1);
E = randn(T,n);
for t = 2:T 
    if t<=ceil(T/2)
    y(t,:) = (A1*y(t-1,:)' + B1*E(t,:)')'   ; 
    z(t) = phi1*E(t,1) + sig_v *randn;
    else
    y(t,:) = (A2*y(t-1,:)' + B2*E(t,:)')' ;
    z(t) = phi1*E(t,1) + sig_v *randn;
    end
    
end
h = 8; alpha = 0.1; 
IRF_1 = zeros(h,n);
IRF_2 = zeros(h,n);
for i = 0:h 
    IRF_1(i+1,:) = (A1^i*B1(:,1))';
    IRF_2(i+1,:) = (A2^i*B2(:,1))'; 
end

X_exo = ones(T,1);
H = T^.5; 
idx = ceil( T/8:(1/8*T):T/8*7 );
idx_frac = idx./T*100;

%% Estimate Bayesian TVP model
tvp_svar_iv_bayes = tvp_svar_iv_bayes( y, z, 2 ,  X_exo(:,2:end), idx, h, [1, 1]);
tvp_svar_iv_bayes.IVSVAR = median(tvp_svar_iv_bayes.store_irf,4);
tvp_svar_iv_bayes.IVSVAR_lb_dm =  quantile(tvp_svar_iv_bayes.store_irf,alpha/2,4);
tvp_svar_iv_bayes.IVSVAR_ub_dm = quantile(tvp_svar_iv_bayes.store_irf,1-alpha/2,4);
%% Estimate Kernel TVP model
[ output_externalIV ] = tvp_svar_iv_kernel(y,  z , p, H , alpha , X_exo, idx , h ); 
%% Estimate Path TVP model
c_max = 50;
[ tvp_svar_iv_alpha_TVP_2 ] = tvp_svar_iv_path(y,  z , p, c_max , alpha , X_exo, idx , h ); 


 
%% Plot results 
color_2 = [128,128,128]./255; 
facealpha = 0.25; 
facealpha2 = 0.25; 
varnames = {'y1','y2'};
xgraph = [(0:h), fliplr((0:h))];  
hfig1 = figure(1); 

a = 1; 
for i = 1:size(idx_frac,2)
    for ii = 1:size(varnames,2)
        subplot(size(idx_frac,2),size(varnames,2),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', output_externalIV.IVSVAR(:,ii,i), 'k:','Linewidth',2);    
        Y_2 = [output_externalIV.IVSVAR_lb_dm(:,ii,i)', fliplr( output_externalIV.IVSVAR_ub_dm(:,ii,i)')  ]; 
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 ); 
        title(varnames(ii),'interpreter','latex')
        if idx(i)<ceil(T/2)
            p3 =  plot((0:h)', IRF_1(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_1(:,ii);
        else
            p3 =  plot((0:h)', IRF_2(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_2(:,ii);
        end  
        if a == 1
            legend([p1,p3],{'TVP (Kernel)','True IRF'},'interpreter','latex')
        end
        xlim([0,h]) 
        if mod(a,size(varnames,2))==1
            ylabel(strcat(num2str(round(idx_frac(i)/100,2)),'T'),'interpreter','latex')
        end   
        a = a + 1; 


    end
end
 
set(gcf,'PaperPositionMode','auto') 
set(hfig1, 'Position', [30 50 1200 800])
hfig1 = tightfig(hfig1) ; 
print(hfig1, 'output/IRFs_simulated_kernel', '-dpdf')
 


%% Create Benchmark Figure 
facealpha = 0.25; 
facealpha2 = 0.25; 
varnames = {'y1','y2'};
xgraph = [(0:h), fliplr((0:h))];   

hfig2 = figure(2);
a = 1; 
for i = 1:size(idx_frac,2)
    for ii = 1:size(varnames,2)
        subplot(size(idx_frac,2),size(varnames,2),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', tvp_svar_iv_bayes.IVSVAR(:,ii,i), 'k:','Linewidth',2);    
        Y_2 = [tvp_svar_iv_bayes.IVSVAR_lb_dm(:,ii,i)', fliplr( tvp_svar_iv_bayes.IVSVAR_ub_dm(:,ii,i)')  ]; 
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 ); 
        title(varnames(ii),'interpreter','latex')
        if idx(i)<ceil(T/2)
            p3 =  plot((0:h)', IRF_1(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_1(:,ii);
        else
            p3 =  plot((0:h)', IRF_2(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_2(:,ii);
        end  
        if a == 1
            legend([p1,p3],{'TVP (Bayes)','True IRF'},'interpreter','latex')
        end
        xlim([0,h]) 
        if mod(a,size(varnames,2))==1
            ylabel(strcat(num2str(round(idx_frac(i)/100,2)),'T'),'interpreter','latex')
        end   
        a = a + 1; 
    end
end
 
set(gcf,'PaperPositionMode','auto') 
set(hfig2, 'Position', [30 50 1200 800])
hfig2 = tightfig(hfig2) ;  
print(hfig2, 'output/IRFs_simulated_bayes', '-dpdf')
 



%% Create Benchmark Figure 
facealpha = 0.25; 
facealpha2 = 0.25; 
varnames = {'y1','y2'};
xgraph = [(0:h), fliplr((0:h))];   


hfig3 = figure(3);
a = 1; 
for i = 1:size(idx_frac,2)
    for ii = 1:size(varnames,2)
        subplot(size(idx_frac,2),size(varnames,2),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', tvp_svar_iv_alpha_TVP_2.IVSVAR(:,ii,i), 'k:','Linewidth',2);    
        Y_2 = [tvp_svar_iv_alpha_TVP_2.IVSVAR_lb_dm(:,ii,i)', fliplr( tvp_svar_iv_alpha_TVP_2.IVSVAR_ub_dm(:,ii,i)')  ]; 
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 ); 
        title(varnames(ii),'interpreter','latex')
        if idx(i)<ceil(T/2)
            p3 =  plot((0:h)', IRF_1(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_1(:,ii);
        else
            p3 =  plot((0:h)', IRF_2(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_2(:,ii);
        end  
        if a == 1
            legend([p1,p3],{'TVP (path)','True IRF'},'interpreter','latex')
        end
        xlim([0,h]) 
        if mod(a,size(varnames,2))==1
            ylabel(strcat(num2str(round(idx_frac(i)/100,2)),'T'),'interpreter','latex')
        end  
        a = a + 1;  
    end
end
 
set(gcf,'PaperPositionMode','auto') 
set(hfig3, 'Position', [30 50 1200 800])
hfig3 = tightfig(hfig3) ;
print(hfig3, 'output/IRFs_simulated_mueller', '-dpdf') 

 

c_max = 5;
[ tvp_svar_iv_alpha_TVP_2 ] = tvp_svar_iv_path(y,  z , p, c_max , alpha , X_exo, idx , h ); 

%% Create Benchmark Figure 
facealpha = 0.25; 
facealpha2 = 0.25; 
varnames = {'y1','y2'};
xgraph = [(0:h), fliplr((0:h))];  
hfig4 = figure(4);
a = 1; 
for i = 1:size(idx_frac,2)
    for ii = 1:size(varnames,2)
        subplot(size(idx_frac,2),size(varnames,2),a)
        hold on; grid on;
        yline(0,'k')
        p1 = plot((0:h)', tvp_svar_iv_alpha_TVP_2.IVSVAR(:,ii,i), 'k:','Linewidth',2);  
        Y_2 = [tvp_svar_iv_alpha_TVP_2.IVSVAR_lb_dm(:,ii,i)', fliplr( tvp_svar_iv_alpha_TVP_2.IVSVAR_ub_dm(:,ii,i)')  ]; 
        bounds2 = fill(xgraph,Y_2,color_2,'LineStyle','none');
        set(bounds2,'FaceColor',color_2,'EdgeColor',color_2,'FaceAlpha',facealpha2 ,'EdgeAlpha',facealpha2 ); 
        title(varnames(ii),'interpreter','latex')
        if idx(i)<ceil(T/2)
            p3 =  plot((0:h)', IRF_1(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_1(:,ii);
        else
            p3 =  plot((0:h)', IRF_2(:,ii), 'k','Linewidth',2);
            IRFi_true = IRF_2(:,ii);
        end  
        if a == 1
            legend([p1,p3],{'TVP (path)','True IRF'},'interpreter','latex')
        end
        xlim([0,h]) 
        if mod(a,size(varnames,2))==1
            ylabel(strcat(num2str(round(idx_frac(i)/100,2)),'T'),'interpreter','latex')
        end   
 
        a = a + 1; 
    end
end
 
set(gcf,'PaperPositionMode','auto') 
set(hfig4, 'Position', [30 50 1200 800])
hfig4 = tightfig(hfig4) ;
print(hfig4, 'output/IRFs_simulated_mueller_cmax_5', '-dpdf') 


