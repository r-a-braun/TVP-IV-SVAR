function [ output ] = tvp_svar_iv_kernel_fast( Y, z, p , H, alpha, X_exo, idx, horizons)
%% Estimates a SVAR-IV's IRFs subject to time-varying parameters.
%  This function (fast) is a short version of tvp_svar_iv_kernel.m used for the simulation study, computing only point estimates and confidence intervals 


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


SN_q = norminv(1-alpha/2,0,1);
Jall = zeros(n*(k-ic)+n+1, n*(k-ic)+n+n*(n+1)/2);
Jall(1:n*(k-ic),1:n*(k-ic)) = eye(n*(k-ic));
Jall(n*(k-ic)+1:n*(k-ic)+n,n*(k-ic)+1:n*(k-ic)+n) = eye(n);

 

%% Loop over each time period specified in the grid
output.waldtest_dof = size(h_sel,1);
output.ftest_dof_1 = [size(h_sel,1),H*n3 - (n3^2)*p - n3];
output.ftest_dof_2 = [size(h_sel,1),H-n3*p-1]; 
a = 1;
for t = idx 
    wt = W(:,a); 
    
    
    
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




    a = a + 1;
end

%% Constant parameter point estimates: 
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
Theta = zeros(n,horizons+1);  
for i = 0:horizons
    Theta(:, i+1) = Phis(:,:,i+1)*B_it; 
end
output.IVSVAR_cp = Theta'; 


end


