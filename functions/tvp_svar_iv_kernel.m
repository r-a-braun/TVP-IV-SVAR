function [ output ] = tvp_svar_iv_kernel( Y, z, p , H, alpha, X_exo, idx, horizons) 
%% Estimates SVAR-IV's IRFs with time-varying parameters (TVP) [Braun, Kapetanios & Marcellino 2025]
% This function estimates Impulse Response Functions (IRFs) for a Structural Vector Autoregressive model 
% with Instrumental Variables (SVAR-IV) using time-varying parameters. The IRFs are estimated at specified 
% time periods (`idx`) and horizons, with uncertainty quantified via Delta Method and Weak IV robust 
% approaches. The estimator is based on the TVP-IV-SVAR methodology,
% assuming a shock of size one standard deviation (absolute IRFs)
% 
% The function also computes F/Wald tests to test invertibility (see text)
% and a series of Portmanteau tests for remaining residual autocorrelations
% It also computes fixed-parameter IV-SVAR estimates 
% Inputs:
%   Y         - T x n matrix of endogenous variables (T: time periods, n: number of variables)
%   z         - T x 1 vector of instrumental variables (IV)
%   p         - Number of lags in the VAR model
%   H         - Bandwidth parameter for the kernel function
%   alpha     - Significance level for constructing confidence intervals (e.g., 0.05 for 95% confidence)
%   X_exo     - T x ic matrix of exogenous variables (intercepts, dummies)
%   idx       - A vector of time indices at which IRFs are to be computed
%   horizons  - Number of periods ahead for which IRFs should be calculated
%
% Outputs:
%   output    - A structure containing the following fields:
%       IVSVAR                - (horizons+1) x n x length(idx) matrix of IRFs
%       IVSVAR_ub_dm          - (horizons+1) x n x length(idx) matrix of upper bounds of IRFs using Delta Method
%       IVSVAR_lb_dm          - (horizons+1) x n x length(idx) matrix of lower bounds of IRFs using Delta Method
%       IVSVAR_ub_osw         - (horizons+1) x n x length(idx) matrix of upper bounds of IRFs using Weak IV robust method
%       IVSVAR_lb_osw         - (horizons+1) x n x length(idx) matrix of lower bounds of IRFs using Weak IV robust method
%       alpha                 - length(idx) x 1 vector of estimated alpha values (for significance testing)
%       alpha_ub              - length(idx) x 1 vector of upper bounds of alpha values
%       alpha_lb              - length(idx) x 1 vector of lower bounds of alpha values
%       Ftest                 - length(idx) x 1 vector of F-test statistics
%       Waldtest              - length(idx) x 1 vector of Wald test statistics
%       Ftest_pval_1          - length(idx) x 1 vector of p-values for F-test (first variant)
%       Ftest_pval_2          - length(idx) x 1 vector of p-values for F-test (second variant)
%       Waldtest_pval         - length(idx) x 1 vector of p-values for Wald test
%       Portmanteau_test_pval - length(idx) x 1 vector of p-values for Portmanteau test
%       Portmanteau_test_testat - length(idx) x 1 vector of test statistics for Portmanteau test
%       Portmanteau_test_dof  - length(idx) x 1 vector of degrees of freedom for Portmanteau test
%       IVSVAR_cp             - (horizons+1) x n x length(idx) matrix of
%       constant parameter (CP) IRFS 
%       IVSVAR_lb_osw_cp       - (horizons+1) x n x length(idx) matrix of lower bounds of IRFs using Weak IV robust method (CP)
%       IVSVAR_ub_osw_cp      - (horizons+1) x n x length(idx) matrix of upper bounds of IRFs using Weak IV robust method (CP)
%       lambda_Wald_CP        - Wald test statistic for Granger causality test using IVs (CP method)
%       lambda_F_CP           - F test statistic for Granger causality test using IVs (CP method)
%       pval_Wald_CP          - p-value for Wald test (CP method)
%       pval_F_CP             - p-value for F test (CP method)
%       pval_F_2_CP           - p-value for second variant of F test (CP method)
%       waldtest_dof          - Degrees of freedom for Wald test
%       ftest_dof_1           - Degrees of freedom for the first variant of the F test
%       ftest_dof_2           - Degrees of freedom for the second variant of the F test 
  
%% Preliminaries
%   Detailed explanation goes here
[T , n] = size(Y);
zp = z(p+1:end);
Lags = lagmatrix(Y,1:p);
ic = size(X_exo,2);
X = [X_exo(p+1:end,:), Lags(p+1:end,:)];
y = Y(p+1:end,:);
k = n*p+ic;
GaussK = @(j,t,H) exp(-0.5*((j-t)/H).^2);
W =  H*GaussK((1:T-p)',idx,H)./sum(GaussK((1:T-p)',idx,H),1);   
 



%% Some preliminaries to use formulas that work with vec(At') opposed to vec(At)
D  = duplication(n); 
L  = elimination( n );  

% Restack test
indicesA = reshape(1:n*k,k,n)';
indicesA = indicesA(:,ic+1:end);
indicesG = (n*k+1:n*k+n)';
indicesS = (n*k+n+1:n*k+n+n*(n+1)/2)';
reidx = [ vec(indicesA); indicesG ; indicesS];

%% Some preliminaries regarding the internal-IV model which the test for invertibility works on
Ya = [z, Y]; 
y3 = Ya(p+1:end,:); 
LagsA = lagmatrix(Ya,1:p);  
X3 = [X_exo(p+1:end,:), LagsA(p+1:end,:)];
k3 = size(X3,2); n3 = n + 1; 
indicesAraw = reshape(1:n3*k3,k3,n3)'; 
indicesA_3 = indicesAraw;
reidxPMvecA = vec(indicesA_3);  
Im = zeros(n3); Im(2:end,1)=1; 
SelM = [zeros(n3,ic),repmat(Im,1,p)]; 
h_sel = eye(length(vec(SelM)));
h_sel(vec(SelM)==0,:) = []; % Select coefficients of vec(A_t') that correspond to lagged z_t  

%% Comopute IRFs at each t in "idx"
output.IVSVAR = zeros(horizons+1,n,length(idx));
output.IVSVAR_ub_dm = zeros(horizons+1,n,length(idx));
output.IVSVAR_lb_dm = zeros(horizons+1,n,length(idx));
output.IVSVAR_ub_osw = zeros(horizons+1,n,length(idx));
output.IVSVAR_lb_osw = zeros(horizons+1,n,length(idx)); 
output.alpha    = zeros(length(idx),1);
output.alpha_ub = zeros(length(idx),1);
output.alpha_lb = zeros(length(idx),1); 
output.Ftest = zeros(length(idx),1);
output.Portmanteau_test_pval = zeros(length(idx),1);
output.Portmanteau_test_testat = zeros(length(idx),1);
output.Portmanteau_test_dof = zeros(length(idx),1);

output.Waldtest = zeros(length(idx),1); 
output.Ftest_pval_1 = zeros(length(idx),1);
output.Ftest_pval_2 = zeros(length(idx),1); 
output.Waldtest_pval = zeros(length(idx),1);

SN_q = norminv(1-alpha/2,0,1);
Jall = zeros(n*(k-ic)+n+1, n*(k-ic)+n+n*(n+1)/2);
Jall(1:n*(k-ic),1:n*(k-ic)) = eye(n*(k-ic));
Jall(n*(k-ic)+1:n*(k-ic)+n,n*(k-ic)+1:n*(k-ic)+n) = eye(n);

 
%% Estimate Constant Model 
A =  (X'*X)\(X'*y); % VAR
Uhat = y - X*A ;
Gamma = 1/T*(Uhat)'*zp;
Sighat = 1/T*(Uhat)'*(Uhat);
sighat = L*vec(Sighat);
iSxx = inv(1./T*X'*X);
eyeiSxx = kron(speye(n),iSxx);
SzxISXX = kron(speye(n), (zp'*X/T)*iSxx); % size(n,n*k)
Xi = zeros(T-p,n*k+n+size(L,1));
for j = 1:T-p
    Xi(j,:) = [vec( X(j,:)'*Uhat(j,:)); vec( Uhat(j,:)'*zp(j) - Gamma); L*vec( Uhat(j,:)'*Uhat(j,:)-Sighat)];
end
St = [ [eyeiSxx,   sparse(n*k,n+size(L,1))]; ...
    [ - SzxISXX, speye(n),   sparse(n,size(L,1))];...
    [sparse(size(L,1),n*k+n), speye(size(L,1))] ];
thetahat    = [ vec(A); Gamma; sighat];
COVhat      = St*(1/T*(Xi)'*(Xi))*St';
RForm.AL    = A(ic+1:end,:)';
RForm.Gamma = thetahat(n*k+1:n*k+n);
RForm.Sigma = reshape(D*thetahat(n*k+n+1:n*k+n+size(L,1)),n,n);
RForm.WHat  = COVhat(reidx,reidx);
RForm.alpha   = sqrt(  RForm.Gamma'*(RForm.Sigma\RForm.Gamma) );
RForm.Gamma_a = [RForm.Gamma; RForm.alpha];
% Jacobian: [beta_t, Gamma_t, alpha_t] < - > [beta_t, gamma_t,sig_t]
Jalpha = 1/2*(RForm.Gamma'*(RForm.Sigma\RForm.Gamma))^(-1/2)*...
    [  2*RForm.Gamma'/RForm.Sigma ,  -kron(RForm.Gamma'/RForm.Sigma,RForm.Gamma'/RForm.Sigma)*D ];
Jall(n*(k-ic)+n+1,n*(k-ic)+1:n*(k-ic)+n+n*(n+1)/2) = Jalpha;
RForm.WHat2 = Jall*RForm.WHat*Jall'; 

%%  Compute IRFs based on Asymptotic Distribution: Delta Method
AA = [RForm.AL;[speye(n*(p-1)),sparse(n*(p-1),n)]];
JJ = [speye(n),sparse(n,n*(p-1))];
Phis = zeros(n,n,horizons+1);
Phis(:,:,1) = eye(n);
Apoweri = AA;
for i = 1:horizons
    Phis(:,:,i+1) = JJ*Apoweri*JJ';
    Apoweri = Apoweri * AA;
end
I1 = eye(n+1); I1(end,:)=[]; I2 = eye(n+1); I2(1:end-1,:)=[];
B_it = (I1*RForm.Gamma_a)./(I2*RForm.Gamma_a);
dBdGamalpha = I1./(I2*RForm.Gamma_a) - [zeros(n),(I1*RForm.Gamma_a)./((I2*RForm.Gamma_a)^2)]; 
Theta = zeros(n,horizons+1);
Theta_var_dm = zeros(n,horizons+1);

for i = 0:horizons
    Theta(:, i+1) = Phis(:,:,i+1)*B_it;
    G_i = zeros(n^2,n^2*p);
    for mm = 0:i-1
        G_i = G_i + kron(JJ*(AA')^(i-1-mm),Phis(:,:,mm+1)) ;
    end
    C_i1 = ( kron(B_it',speye(n))*G_i );
    Cbar_i =  Phis(:,:,i+1)*dBdGamalpha ;
    GradAG2 = [C_i1, Cbar_i];
    COVthetat = GradAG2*RForm.WHat2*GradAG2';
    Theta_var_dm(:,i+1) =  diag( COVthetat );

    %%  Compute IRFs based on Asymptotic Distribution: Montiel Olea, Stock & Watson
    C_i1 = ( kron((I1*RForm.Gamma_a)',speye(n))*G_i );
    dLam = [C_i1, Phis(:,:,i+1)*I1  ];
    GradAG_AR2 = [dLam;[zeros(1,size(GradAG2,2)-1),1]];
    COVthetat_AR2 = GradAG_AR2*RForm.WHat2*GradAG_AR2';
    e = eye(n);  ep1 = eye(n+1);
    critval = norminv(1-alpha/2,0,1)^2;
    CIs2 = zeros(2,n);
    for ii = 1:n
        XYf =  [e(ii,:)*Phis(:,:,i+1)*(I1*RForm.Gamma_a); I2*RForm.Gamma_a];
        What = ep1([ii,n+1],:)*COVthetat_AR2*ep1([ii,n+1],:)'/T;
        f0 = XYf(1)^2-critval*What(1,1);
        f1 = XYf(1)*XYf(2) - critval*What(1,2);
        f2 = XYf(2)^2 - critval*What(2,2);
        Dar = f1^2 - f0*f2;
        r1 = (f1-Dar^(1/2))/f2;
        r2 = (f1+Dar^(1/2))/f2;
        if and(Dar>=0,f2>=0)
            CIs2(:,ii) = [r1,r2];
        elseif and(Dar<0,f2<0)
            CIs2(:,ii) = [-inf,inf];
        elseif and(Dar>=0,f2<0)
            CIs2(:,ii) =  [-inf,inf];
        end
    end
    output.IVSVAR_lb_osw_cp(i+1,: ) = CIs2(1,:)';
    output.IVSVAR_ub_osw_cp(i+1,: ) = CIs2(2,:)';
end
output.IVSVAR_cp = Theta';  
Im = zeros(n3); Im(2:end,1)=1; 
SelM = [zeros(n3,ic),repmat(Im,1,p)]; 
h_sel = eye(length(vec(SelM)));
h_sel(vec(SelM)==0,:) = []; % Select coefficients of vec(A_t') that correspond to lagged z_t   
[output.lambda_Wald_CP, output.lambda_F_CP, output.pval_Wald_CP, output.pval_F_CP , output.pval_F_2_CP] = Granger_Test(Ya, p, X_exo, h_sel );

hQ = p*2;
output.Portmanteau_test_h = hQ;
Fq = [];
for i = 1:hQ
    I1 = [zeros(i,T-p-i); eye(T-p-i)];
    I2 = [eye(T-p-i); zeros(i,T-p-i)];
    Fq = [Fq,  I1*I2' ];
end

%% Loop over each time period specified in the grid
output.waldtest_dof = size(h_sel,1);
output.ftest_dof_1 = [size(h_sel,1),H*n3 - (n3^2)*p - n3];
output.ftest_dof_2 = [size(h_sel,1),H-n3*p-1]; 
a = 1;
for t = idx 
    wt = W(:,a);  
    %% Compute the Test for Invertibility 
    wX3 = wt.*X3;
    At3 =  (wX3'*X3)\(wX3'*y3); % VAR-augmented 
    Uhat3 = y3 - X3*At3; 
    Sighat3 = 1/H*(Uhat3.*wt)'*(Uhat3);   
    iPix = inv(1./H*wX3'*X3); Piww =1./H*wX3'*wX3;  
    Asym_COVhatbeta = kron(Sighat3,iPix*Piww*iPix);     
    sig_b = Asym_COVhatbeta(reidxPMvecA,reidxPMvecA)/H;   
    vbhat = vec(At3');
    wald_test = (h_sel*vbhat)'*((h_sel*sig_b*h_sel')\(h_sel*vbhat));
    F_test = wald_test/size(h_sel,1);
    p_wtest = 1-chi2cdf(wald_test,size(h_sel,1));
    p_ftest_1 = 1-fcdf(F_test,size(h_sel, 1), H*n3 - (n3^2)*p - n3); 
    p_ftest_2 = 1-fcdf(F_test,size(h_sel, 1), H-n3*p-1);   
    output.Waldtest(a) = wald_test;
    output.Waldtest_pval(a) = p_wtest;  
    output.Ftest(a) = F_test;
    output.Ftest_pval_1(a) = p_ftest_1;
    output.Ftest_pval_2(a) = p_ftest_2; 
 
    %% SVAR-IV 
    wX = wt.*X;
    At =  (wX'*X)\(wX'*y); % VAR
    Uhat = y - X*At;
    Gamma = 1/H*(Uhat.*wt)'*zp;
    Sighat = 1/H*(Uhat.*wt)'*(Uhat);
    sighat = L*vec(Sighat);
    iSxx = inv(1./H*wX'*X);
    eyeiSxx = kron(speye(n),iSxx);
    SzxISXX = kron(speye(n), (zp'*wX/H)*iSxx); % size(n,n*k)
    Xi = zeros(T-p,n*k+n+size(L,1));
    for j = 1:T-p
        Xi(j,:) = [vec( X(j,:)'*Uhat(j,:)); vec( Uhat(j,:)'*zp(j) - Gamma); L*vec( Uhat(j,:)'*Uhat(j,:)-Sighat)];
    end
    St = [ [eyeiSxx,   sparse(n*k,n+size(L,1))]; ...
        [ - SzxISXX, speye(n),   sparse(n,size(L,1))];...
        [sparse(size(L,1),n*k+n), speye(size(L,1))] ];
    thetahat    = [ vec(At); Gamma; sighat];
    COVhat      = St*(1/H*(Xi.*wt)'*(Xi.*wt))*St'; 
    RForm.AL    = At(ic+1:end,:)';
    RForm.Gamma = thetahat(n*k+1:n*k+n);
    RForm.Sigma = reshape(D*thetahat(n*k+n+1:n*k+n+size(L,1)),n,n);
    RForm.WHat  = COVhat(reidx,reidx);    
    RForm.alpha   = sqrt(  RForm.Gamma'*(RForm.Sigma\RForm.Gamma) );
    RForm.Gamma_a = [RForm.Gamma; RForm.alpha];
    % Jacobian: [beta_t, Gamma_t, alpha_t] < - > [beta_t, gamma_t,sig_t]
    Jalpha = 1/2*(RForm.Gamma'*(RForm.Sigma\RForm.Gamma))^(-1/2)*...
        [  2*RForm.Gamma'/RForm.Sigma ,  -kron(RForm.Gamma'/RForm.Sigma,RForm.Gamma'/RForm.Sigma)*D ];   
    Jall(n*(k-ic)+n+1,n*(k-ic)+1:n*(k-ic)+n+n*(n+1)/2) = Jalpha;   
    RForm.WHat2 = Jall*RForm.WHat*Jall';   
    alpha_se = sqrt(RForm.WHat2(end,end)/H);
    output.alpha(a)    = RForm.alpha;
    output.alpha_ub(a) = RForm.alpha + alpha_se*SN_q;
    output.alpha_lb(a) = RForm.alpha - alpha_se*SN_q;   
      
    
    
    
    
    %%  Compute IRFs based on Asymptotic Distribution: Delta Method
    AA = [RForm.AL;[speye(n*(p-1)),sparse(n*(p-1),n)]];
    JJ = [speye(n),sparse(n,n*(p-1))];
    Phis = zeros(n,n,horizons+1);
    Phis(:,:,1) = eye(n);
    Apoweri = AA;
    for i = 1:horizons
        Phis(:,:,i+1) = JJ*Apoweri*JJ';
        Apoweri = Apoweri * AA;
    end
    I1 = eye(n+1); I1(end,:)=[]; I2 = eye(n+1); I2(1:end-1,:)=[];
    B_it = (I1*RForm.Gamma_a)./(I2*RForm.Gamma_a);
    dBdGamalpha = I1./(I2*RForm.Gamma_a) - [zeros(n),(I1*RForm.Gamma_a)./((I2*RForm.Gamma_a)^2)];  
    
    Theta = zeros(n,horizons+1);
    Theta_var_dm = zeros(n,horizons+1); 
    
    for i = 0:horizons
        Theta(:, i+1) = Phis(:,:,i+1)*B_it;
        G_i = zeros(n^2,n^2*p);
        for mm = 0:i-1
            G_i = G_i + kron(JJ*(AA')^(i-1-mm),Phis(:,:,mm+1)) ;
        end
        C_i1 = ( kron(B_it',speye(n))*G_i );
        Cbar_i =  Phis(:,:,i+1)*dBdGamalpha ;
        GradAG2 = [C_i1, Cbar_i];
        COVthetat = GradAG2*RForm.WHat2*GradAG2'; 
        Theta_var_dm(:,i+1) =  diag( COVthetat );
        
        %%  Compute IRFs based on Asymptotic Distribution: Montiel Olea, Stock & Watson 
        C_i1 = ( kron((I1*RForm.Gamma_a)',speye(n))*G_i );
        dLam = [C_i1, Phis(:,:,i+1)*I1  ]; 
        GradAG_AR2 = [dLam;[zeros(1,size(GradAG2,2)-1),1]];
        COVthetat_AR2 = GradAG_AR2*RForm.WHat2*GradAG_AR2'; 
        e = eye(n);  ep1 = eye(n+1);
        critval = norminv(1-alpha/2,0,1)^2; 
        CIs2 = zeros(2,n);
        for ii = 1:n  
            XYf =  [e(ii,:)*Phis(:,:,i+1)*(I1*RForm.Gamma_a); I2*RForm.Gamma_a];
            What = ep1([ii,n+1],:)*COVthetat_AR2*ep1([ii,n+1],:)'/H;
            f0 = XYf(1)^2-critval*What(1,1);
            f1 = XYf(1)*XYf(2) - critval*What(1,2);
            f2 = XYf(2)^2 - critval*What(2,2);
            Dar = f1^2 - f0*f2;
            r1 = (f1-Dar^(1/2))/f2;
            r2 = (f1+Dar^(1/2))/f2;
            if and(Dar>=0,f2>=0)
                CIs2(:,ii) = [r1,r2];
            elseif and(Dar<0,f2<0)
                CIs2(:,ii) = [-inf,inf];
            elseif and(Dar>=0,f2<0)
                CIs2(:,ii) =  [-inf,inf];
            end
        end 
        output.IVSVAR_lb_osw(i+1,:,a) = CIs2(1,:)';
        output.IVSVAR_ub_osw(i+1,:,a) = CIs2(2,:)';  
    end 
    output.IVSVAR(:,:,a) = Theta';
    IRF_sd2 = sqrt( Theta_var_dm ./H );
    output.IVSVAR_ub_dm(:,:,a) = (Theta + IRF_sd2.*SN_q)';
    output.IVSVAR_lb_dm(:,:,a) = (Theta - IRF_sd2.*SN_q)';




    %% Residual autocorrelation tests:
    Uhatp = Uhat';
    Uhatpw = (Uhatp.*wt');  
    C_h = 1/H*(Uhatp*Fq)*kron(eye(hQ),Uhatpw');
    iC_0 = inv(Uhatp*Uhatpw'/H);
    Qh = H*vec(C_h)'*kron(eye(hQ),kron(iC_0,iC_0))*vec(C_h)  ;
    dofQ = n^2*(hQ-p); 
    p_Qtest = 1-chi2cdf(Qh,dofQ) ; 
    output.Portmanteau_test_pval(a) = p_Qtest;
    output.Portmanteau_test_testat(a) = Qh;
    output.Portmanteau_test_dof(a) = dofQ;


    a = a + 1;
end


end


