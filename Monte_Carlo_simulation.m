%% This code replicates the simulation results
clear; clc;
rng(100)
addpath('functions','data')
p = 3; % Lag length VAR(p)
factorT = 1; % Can be used to scale up T as described in the online appendix

%% Step 1: Estimate the DGP (see Montiel Olea et al [2021])
data_raw = table2array(readtable('data_kilian_09.txt'));
data = data_raw(1:380,:); %  1973:M1–2004:M9
Htrue = 100;
Y0 = data(1:p,:);
shortY = data(p+1:end,:);
[Tdata, n] = size(shortY);
X = lagmatrix(data,1:p);
X = [ones(Tdata,1),X(p+1:end,:)];
Bols = (X'*X)\(X'*shortY);
U = shortY-X*Bols;
Sigma = U'*U/Tdata;
eaux  = [1;1;-1];    % B1
B1 = eaux*((eaux'*(Sigma)^(-1)*eaux)^(-.5)); %
q = inv(chol(Sigma,'lower'))*B1;
Q = [q,null(q')]; % Rotation matrix

GaussK = @(j,t,Hk) exp(-0.5*((j-t)/Hk).^2);
Wdata =  Htrue*GaussK((1:Tdata)',(1:Tdata),Htrue)./sum(GaussK((1:Tdata)',(1:Tdata),Htrue),1);
for t = 1:Tdata
    wt = Wdata(:,t);
    wX = wt.*X;
    At =  (wX'*X)\(wX'*shortY);
    Uhat = shortY - X*At;
    Sigt = 1/Htrue*(Uhat.*wt)'*(Uhat);
    Bt = chol(Sigt,'lower')*Q;
    vparamTraw(:,t) = [vec(At); vec(Bt)] ;
end

%% DO SIMU
T = floor(factorT*Tdata); idx_inter = (1:T)./factorT; idx_inter(idx_inter<1)=1;
vparamT = zeros(size(vparamTraw,1),length(idx_inter));
for i = 1:size(vparamTraw,1)
    vparamT(i,:) = interp1((1:Tdata),vparamTraw(i,:),idx_inter);
end

% Compute the True Impulse Response Functions:
h = 20;
t_stand = ceil(T/2);
Btstand = reshape(vparamT((n*(n*p+1))+1:end,t_stand),n,n);
IRFtrue_bar = zeros(h+1,n,T);
IRFtrue_sig = zeros(h+1,n,T);
for t = 1:T
    vAvB = vparamT(:,t);
    At = reshape(vAvB(1:(n*(n*p+1))),n*p+1,n);
    Bt = reshape(vAvB((n*(n*p+1))+1:end),n,n);
    [ AA, J , nu, Afull] = Companion( At' , 1);
    Btbar = Bt(:,1)./Btstand(1,1);
    irftsig(1,:) = Bt(:,1);
    irftbar(1,:) = Btbar;
    for j = 1:h
        irftsig(j+1,:) = J*AA^j*J'*Bt(:,1);
        irftbar(j+1,:) = J*AA^j*J'*Btbar;
    end
    IRFtrue_bar(:,:,t) = irftbar;
    IRFtrue_sig(:,:,t) = irftsig;
end
yvars = {'{dprod}_t','{rea}_t','{rpo}_t'};
% Figure 1: True impulse response functions
hfig = figure(1); 
time_index = (1:T)'; idx_table =  ceil(linspace(1,T,20)') ;
T_figure_all = table(time_index(idx_table));
T_figure_all.Properties.VariableNames = {'Time'};
for i = 1:size(data,2)
    subplot(1,size(data,2),i)
    hold on; grid on;
    p0 = plot(1:T,squeeze(IRFtrue_sig(1,i,:)),'k','linewidth',2) ;
    p3 = plot(1:T,squeeze(IRFtrue_sig(1+10,i,:)),'k:','linewidth',2);
    p5 = plot(1:T,squeeze(IRFtrue_sig(1+20,i,:)),'k--','linewidth',2);
    xlim([1,T])
    xlabel('Time'); ylabel('IRF')
    if i==1
        legend([p0 p3 p5],{'h=0','h=10','h=20'},'Interpreter','Latex')
    end
    title(strcat('T= ',num2str(T), ',\ $\varepsilon_{1t} \rightarrow \, ',yvars{i},'$'),'interpreter','latex')
end
set(hfig,'PaperPositionMode','auto')
set(hfig, 'Position', [30 50 800 400])
hfig = tightfig(hfig);
print(hfig, 'output/simulation_true_irfs' , '-dpdf')

%% Step 2: Simulate from the DGP and estimate confidence intervals.
alpha = 0.05; % Confidence intervals
for estimate_H = [0,1]
    for weakIV = [0, 1]
        if weakIV == 0
            phiz = 0.8600;
            sigz = 0.0638;
        elseif weakIV == 1
            phiz = 0.4818;
            sigz = 0.7152;
        end
        nrep = 4999;   n = size(data,2);
        idx = [ceil(T/2), ceil(T*3/4)]; % Compute Confidence Intervals at 2 points of time
        npoints = length(idx);

        IRFtrue_bar_idx = IRFtrue_bar(:,:,idx);
        IRFtrue_sig_idx = IRFtrue_sig(:,:,idx);

        % Relative IRFs based on internal instrument VAR
        lambda_tilde_lb_dm = zeros(h+1,n,npoints,nrep);
        lambda_tilde_ub_dm = zeros(h+1,n,npoints,nrep);
        lambda_tilde_lb_ar = zeros(h+1,n,npoints,nrep);
        lambda_tilde_ub_ar = zeros(h+1,n,npoints,nrep);

        % Absolute IRFs based on IV-SVAR
        lambda_lb_dm = zeros(h+1,n,npoints,nrep);
        lambda_ub_dm = zeros(h+1,n,npoints,nrep);
        lambda_lb_ar = zeros(h+1,n,npoints,nrep);
        lambda_ub_ar = zeros(h+1,n,npoints,nrep);


        ic = 1; idx_shock = 1 ;
        for i = 1:nrep
            clc; disp(i/nrep)
            % Simulate the TVP-VAR using the estimated Law of Motions
            nburn = 100;
            Ysim = zeros(T+p+nburn,n); usim = zeros(T+nburn,n);
            esim = zeros(T+nburn,n); Ysim(1:p,:) = Y0(1:p,:);
            tsim = 1;
            vparamTb = [repmat(vparamT(:,1),1,nburn),vparamT];
            for t = p+1:T+p+nburn
                vAvB = vparamTb(:,t-p);
                At = reshape(vAvB(1:(n*(n*p+1))),n*p+1,n);
                Bt = reshape(vAvB((n*(n*p+1))+1:end),n,n);
                [ AA, J , nu, Afull] = Companion(At' , 1);
                Ylag = vec(Ysim(t-1:-1:t-p,:)');
                esim(tsim,:) = randn(n,1);
                usim(tsim,:) = (Bt*esim(tsim,:)')';
                Ysim(t,:) = (J*nu)' +  (J*(AA*Ylag)+usim(tsim,:)')' ;
                tsim = tsim + 1;
            end
            zsim = [zeros(p,1); phiz.*esim(:,1) + sigz.*randn(T+nburn,1)];
            Ysim = Ysim(nburn+1:end,:); zsim = zsim(nburn+1:end,:);
            X_exo = ones(T+p,1);
            if estimate_H == 1
               Hgrid = T.^(0.5:0.01:.8);  
               [Hest, CVs ]  = choose_H_CF([zsim, Ysim],p,1,Hgrid,ceil(factorT));   
            else
                Hest = floor(Htrue*sqrt(factorT)); 
            end
            [ output_lam_ ] = tvp_svar_iv_kernel( Ysim, zsim, p , Hest, alpha, X_exo, idx, h);
            [ output_lam_til ] = tvp_svar_internal_kernel( Ysim, zsim, p , Hest, alpha, X_exo, idx, h, idx_shock, t_stand);

            %% Save results
            lambda_lb_dm(:,:,:,i) = output_lam_.IVSVAR_lb_dm;
            lambda_ub_dm(:,:,:,i) = output_lam_.IVSVAR_ub_dm;
            lambda_lb_ar(:,:,:,i) = output_lam_.IVSVAR_lb_osw;
            lambda_ub_ar(:,:,:,i) = output_lam_.IVSVAR_ub_osw;

            lambda_tilde_lb_dm(:,:,:,i) = output_lam_til.VARpm_lb_dm;
            lambda_tilde_ub_dm(:,:,:,i) = output_lam_til.VARpm_ub_dm;
            lambda_tilde_lb_ar(:,:,:,i) = output_lam_til.VARpm_lb_msw;
            lambda_tilde_ub_ar(:,:,:,i) = output_lam_til.VARpm_ub_msw; 

        end


        %  Compute Empirical Coverage
        CI_lambda_DM = zeros(h+1,n,npoints,nrep);
        CI_lambda_AR = zeros(h+1,n,npoints,nrep);
        CI_lambda_til_DM = zeros(h+1,n,npoints,nrep);
        CI_lambda_til_AR = zeros(h+1,n,npoints,nrep);
        for i = 1:nrep
            CI_lambda_DM(:,:,:,i) = and(lambda_lb_dm(:,:,:,i)<=IRFtrue_sig_idx,IRFtrue_sig_idx<=lambda_ub_dm(:,:,:,i));
            CI_lambda_AR(:,:,:,i) = and(lambda_lb_ar(:,:,:,i)<=IRFtrue_sig_idx,IRFtrue_sig_idx<=lambda_ub_ar(:,:,:,i));
            CI_lambda_til_DM(:,:,:,i) = and(lambda_tilde_lb_dm(:,:,:,i)<=IRFtrue_bar_idx,IRFtrue_bar_idx<=lambda_tilde_ub_dm(:,:,:,i));
            CI_lambda_til_AR(:,:,:,i) = and(lambda_tilde_lb_ar(:,:,:,i)<=IRFtrue_bar_idx,IRFtrue_bar_idx<=lambda_tilde_ub_ar(:,:,:,i));
        end 
        EmpCovRateLam_DM = mean(CI_lambda_DM,4);
        EmpCovRateLam_AR = mean(CI_lambda_AR,4);
        EmpCovRateLamTil_AR = mean(CI_lambda_til_AR,4);
        EmpCovRateLamTil_DM = mean(CI_lambda_til_DM,4);
 


        %% Step 3: Plot the output
        xvals = (0:h)'; xlinecov = ones(h+1,1).*(1-alpha); idx_plot = (0:5:20)+1;
        hfig0 = figure(200);
        a = 1;
        for ii = 1:npoints
            for i = 1:n
                subplot(npoints,2*n,a)
                hold on; grid on;
                p1 = plot(xvals(idx_plot),EmpCovRateLam_DM(idx_plot,i,ii),'k--diamond','linewidth',1.2) ;
                p2 = plot(xvals(idx_plot),EmpCovRateLam_AR(idx_plot,i,ii),'k--*','linewidth',1.2) ;
                p3 = plot(xvals(idx_plot),xlinecov(idx_plot),'k');

                xlim([0,h]);
                ylim([.7,1])
                if ii == npoints
                    xlabel('horizons','interpreter','latex');
                end
                if or(a == 1,a==npoints*n+1)
                    ylabel('coverage','interpreter','latex');
                end
                if a<=3
                    title({'absolute IRF',strcat('$\lambda_{h,',num2str(i),',',num2str(idx(ii)),'}$')},'interpreter','latex');
                else
                    title(strcat('$\lambda_{h,',num2str(i),',',num2str(idx(ii)),'}$'),'interpreter','latex');
                end
                if a == 9
                    legend([p3 p1 p2],{'nom. Coverage','$CS^{DM}(95\%)$','$CS^{AR}(95\%)$'},'interpreter','latex','location','southeast')
                end
                xticks([0,10,20])
                a = a + 1;
            end
            for i = 1:n
                subplot(npoints,2*n,a)
                hold on; grid on;
                p1 = plot(xvals(idx_plot),EmpCovRateLamTil_DM(idx_plot,i,ii),'k--diamond','linewidth',1.2);
                p2 = plot(xvals(idx_plot),EmpCovRateLamTil_AR(idx_plot,i,ii),'k--*','linewidth',1.2);
                p3 = plot(xvals(idx_plot),xlinecov(idx_plot),'k');
                xlim([0,h]);
                ylim([.7,1])
                if ii == npoints
                    xlabel('horizons','interpreter','latex');
                end
                if a<=9
                    title({'relative IRF',strcat('$\tilde{\lambda}_{h,',num2str(i),',',num2str(idx(ii)),'}$')},'interpreter','latex');
                else
                    title(strcat('$\tilde{\lambda}_{h,',num2str(i),',',num2str(idx(ii)),'}$'),'interpreter','latex');
                end
                a = a + 1;
                xticks([0,10,20])
            end
            if and(ii == npoints,i == n)
                legend([p3 p1 p2],{'nom. Coverage','$CS^{DM}(95\%)$','$CS^{AR}(95\%)$'},'interpreter','latex',...
                    'location','southeast')
            end
        end
        set(hfig0,'PaperPositionMode','auto')
        set(hfig0, 'Position', [30 50 700 300])
         
        hfig0 = tightfig(hfig0);
        if estimate_H == 0
            print(hfig0,strcat('output/fig_sim_Hknown_',num2str(Htrue),'_T_',num2str(T),'_weakIV_',num2str(weakIV)), '-dpdf')
        elseif estimate_H == 1
            print(hfig0,strcat('output/fig_sim_Hest_',num2str(Htrue),'_T_',num2str(T),'_weakIV_',num2str(weakIV)), '-dpdf')
        end

        close all


    end


end






