function output = tvp_svar_iv_path( Y, ziv, p , c_max, alpha, X_exo, idx, horizons)
% Function to compute the Mueller Petalas Path estimator


%% Preliminaries
%   Detailed explanation goes here
[TpP , n] = size(Y);
T = TpP-p;
zp = ziv(p+1:end);
Lags = lagmatrix(Y,1:p);
ic = size(X_exo,2);
X = [X_exo(p+1:end,:), Lags(p+1:end,:)]';
y = Y(p+1:end,:);
Y0 = Y(p+1:end,:)';   
k = n*p+ic;   
Ahat = (Y0*X')/(X*X');
Uhat = Y0 - Ahat*X; 
SN_q = norminv(1-alpha/2,0,1);

%% Some preliminaries to use formulas that work with vec(At') opposed to vec(At)
D  = duplication(n) ; 
L  = elimination( n+1 );    
%% Compute Gradient and Hessian 
Sighat = Uhat*Uhat'/T;
iSighat = inv(Sighat);

Uhat_a = [ zp';Uhat];
Sighat_a = Uhat_a*Uhat_a'/T;
iSighat_a = inv(Sighat_a);
%Uhat*zp/T
indicesSigma = tril(ones(n+1));
indicesSigma(:,1) = 0; 

sighata = L*vec(Sighat_a);
ahat = vec(Ahat); 
thetahat = [ahat;sighata]; 
grad_a = zeros(length(ahat),T); 
grad_sig = zeros(length(sighata),T);  

for t = 1:T
    % Gradient with respect to "a"
    grad_a(:,t) = kron(X(:,t),iSighat)*Y0(:,t) - kron(X(:,t)*X(:,t)',iSighat)*ahat;
    % Gradient with respect to "vech(Sigma)"
    grad_sig(:,t) = L*vec(-1/2*iSighat_a + 1/2*iSighat_a*Uhat_a(:,t)*Uhat_a(:,t)'*iSighat_a);
end
H_a = - kron(X*X',iSighat);
H_s = L*(T/2*kron(iSighat_a,iSighat_a)...
    - 1/2*kron(iSighat_a,iSighat_a*Uhat_a*Uhat_a'*iSighat_a)...
    - 1/2*kron(iSighat_a*Uhat_a*Uhat_a'*iSighat_a,iSighat_a))*L';


sgrad = [ grad_a ; grad_sig];
H = - blkdiag(H_a, H_s) / T; % Hessian is asymptotically block diagonal since U'*X/T -> 0
V =  sgrad*sgrad'/T; % Covariance of Score 
 
S =  (H\V)/H ; % Sandwich Estimator 

%% Estimate the constant parameter IRFs
S_gamma = eye((n+2)*(n+1)/2);
S_gamma([1,n+2:(n+2)*(n+1)/2],:)=[];

S_sigma = eye((n+2)*(n+1)/2);
S_sigma(L*vec(indicesSigma)~=1,:)=[];
indicesA = reshape(1:n*k,n,k) ;
indicesA = indicesA(:,ic+1:end);
reidx = [ vec(indicesA); n*k+find(sum(S_gamma,1))' ; n*k+find(sum(S_sigma,1))'];
idx_A = vec(indicesA);
idx_Gamma = n*k+find(sum(S_gamma,1))' ;
idx_Sigma = n*k+find(sum(S_sigma,1))';

RForm.AL    = Ahat(:,ic+1:end);
RForm.Gamma = S_gamma*sighata;
RForm.Sigma = reshape(D*S_sigma*sighata,n,n);
RForm.WHat  = S(reidx,reidx);

% Jacobian: [beta_t, Gamma_t, alpha_t] < - > [beta_t, gamma_t,sig_t]
Jalpha = 1/2*(RForm.Gamma'*(RForm.Sigma\RForm.Gamma))^(-1/2)*...
    [  2*RForm.Gamma'/RForm.Sigma ,  -kron(RForm.Gamma'/RForm.Sigma,RForm.Gamma'/RForm.Sigma)*D ];
Jall = zeros(n*(k-ic)+n+1, n*(k-ic)+n+n*(n+1)/2);
Jall(1:n*(k-ic),1:n*(k-ic)) = eye(n*(k-ic));
Jall(n*(k-ic)+1:n*(k-ic)+n,n*(k-ic)+1:n*(k-ic)+n) = eye(n);
Jall(n*(k-ic)+n+1,n*(k-ic)+1:n*(k-ic)+n+n*(n+1)/2) = Jalpha;
RForm.WHat2 = Jall*RForm.WHat*Jall'; 
RForm.alpha   = sqrt(  RForm.Gamma'*(RForm.Sigma\RForm.Gamma) );
RForm.Gamma_a = [RForm.Gamma; RForm.alpha];


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
%  [eye(n)./RForm.Gamma_a(end), - RForm.Gamma_a(1:end-1)./RForm.Gamma_a(end)^2]
% [zeros(n),(I1*RForm.Gamma_a)./((I2*RForm.Gamma_a)^2)]

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

%% Start the Mueller and Petalas algorithm for the TVP reduced form estimates
% and covariance matrix
C = linspace(.1,c_max,10); n_G = size(C,2);
xt = H\sgrad;
dxt = diff(xt')';
ytil = H*(V\sgrad);
p_ = length(thetahat); % # of parameters
wtil = NaN(n_G,1); logwtil =  NaN(n_G,1); logwtil2 =  NaN(n_G,1);
betas_all = zeros(p_,T,n_G);
qLL_ci = NaN(n_G,1);
for i = 1:n_G

 

    c_i = C(i); 
    r_i = 1 - c_i/T;
    z = zeros(size(xt,1),T);
    z(:,1) = xt(:,1);
    for t = 2:T
        z(:,t) = r_i*z(:,t-1) + dxt(:,t-1); %xt(:,t)-xt(:,t-1)
    end 
    
    Xrhs2 = (r_i.^( (1:T)-1))';
    beta_reg2 = (Xrhs2'*Xrhs2)\(Xrhs2'*z');
    ztil = (z'- Xrhs2*beta_reg2)';

    % (c) Create zbar
    zbar = zeros(size(ztil,1),T);
    zbar(:,T) = ztil(:,T);
    for t = (T-1):-1:1
        zbar(:,t) = r_i*zbar(:,t+1) + ztil(:,t)- ztil(:,t+1);
    end
    % (d) create betahat_i
    betas_all(:,:,i) = thetahat + xt - r_i*zbar;
 

    % (E) compute mixture probability
    qLL_ci2 = 0;
    for t= 1:T
        qLL_ci2 = qLL_ci2 + (r_i*zbar(:,t) - xt(:,t))'*ytil(:,t);
        
    end

    qLL_ci(i) = sum(diag( ((r_i*zbar - xt)'*ytil) ) ); 
    if i == 1
        wtil(i) = 1; 
    else
        wtil(i) = sqrt( T*(1-r_i^2)*r_i^(T-1)/(1-r_i^(2*T)) )*exp(-1/2*qLL_ci(i)); 
    end
    if i == 1
        logwtil(i) = 1; 
    else
        logwtil(i) = 0.5*(log(T)+log(1-r_i^2)+(T-1)*log(r_i) -  log(1-r_i^(2*T))) - 1/2*qLL_ci(i) ;

    end
    

    %disp('stop')
end
w = (wtil./sum(wtil))';  
maxll = max(logwtil);
w = exp(logwtil - maxll - log(sum( exp(logwtil-maxll) ))     )'  ;  

w3d = reshape(w, 1, 1, []);
betahat = sum(betas_all.*w3d, 3);
 
%% Compute Covariance estimate at a few time points specified in idx 
ktc = @(t,c) ( c*(1+exp(2*c)+exp(2*c*t/T)+exp(2*c*(1-t/T))) )...
    /(2*exp(2*c)-2);



%% Loop over each time period specified in the grid
a = 1;
for t = idx 
    % Compute covariance matrix:
    Omega_t = zeros(size(betahat,1),size(betahat,1));
    for j = 1:n_G
        c_i = C(j);
        if c_i == 0
            kt_c = 1;
        else
            kt_c = ktc(t,c_i);
        end
        Omega_t = Omega_t + w(j)*( S/T*kt_c + ...
            (betas_all(:,t,j)-betahat(:,t))*(betas_all(:,t,j)-betahat(:,t))' ); 
    end
    vecASiga_t = betahat(:,t); 
    RForm.AL    = reshape(vecASiga_t(idx_A),n,k-ic); 
    RForm.Gamma = vecASiga_t(idx_Gamma);
    RForm.Sigma = reshape(D*vecASiga_t(idx_Sigma),n,n);
    RForm.WHat  = Omega_t(reidx,reidx);
    % Jacobian: [beta_t, Gamma_t, alpha_t] < - > [beta_t, gamma_t,sig_t]
    Jalpha = 1/2*(RForm.Gamma'*(RForm.Sigma\RForm.Gamma))^(-1/2)*...
        [  2*RForm.Gamma'/RForm.Sigma ,  -kron(RForm.Gamma'/RForm.Sigma,RForm.Gamma'/RForm.Sigma)*D ]; 
    Jall(n*(k-ic)+n+1,n*(k-ic)+1:n*(k-ic)+n+n*(n+1)/2) = Jalpha;
    RForm.WHat2 = Jall*RForm.WHat*Jall';
    RForm.alpha   = sqrt(  RForm.Gamma'*(RForm.Sigma\RForm.Gamma) );
    RForm.Gamma_a = [RForm.Gamma; RForm.alpha];

    RForm.WHat2 = Jall*RForm.WHat*Jall';
    alpha_se = sqrt(RForm.WHat2(end,end));
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
  %  [eye(n)./RForm.Gamma_a(end), - RForm.Gamma_a(1:end-1)./RForm.Gamma_a(end)^2] 
   % [zeros(n),(I1*RForm.Gamma_a)./((I2*RForm.Gamma_a)^2)] 
    
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
            What = ep1([ii,n+1],:)*COVthetat_AR2*ep1([ii,n+1],:)';
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
    IRF_sd2 = sqrt( Theta_var_dm  );
    output.IVSVAR_ub_dm(:,:,a) = (Theta + IRF_sd2.*SN_q)';
    output.IVSVAR_lb_dm(:,:,a) = (Theta - IRF_sd2.*SN_q)';
    a = a + 1;
end

end
 