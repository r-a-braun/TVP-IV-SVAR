
function [ output ] = tvp_svar_internal_kernel_fast( Y, z, p , H, alpha, X_exo, idx, horizons, idx_shock, t_standard)
%% Estimates a SVAR-IV's IRFs subject to time-varying parameters (relative IRFs based on internal instrument)
%  This function (fast) is a short version of tvp_svar_internal_kernel_fast.m used for the simulation study, computing only point estimates and confidence intervals 


%% Preliminaries
%   Detailed explanation goes here
[T , n] = size(Y);
Ya = [z, Y];  
LagsA = lagmatrix(Ya,1:p);   
ic = size(X_exo,2); 
X3 = [X_exo(p+1:end,:), LagsA(p+1:end,:)]; 
y3 = Ya(p+1:end,:);  
k3 = size(X3,2);
GaussK = @(j,t,H) exp(-0.5*((j-t)/H).^2); 
W =  H*GaussK((1:T-p)',idx,H)./sum(GaussK((1:T-p)',idx,H),1);   



%% Some preliminaries to use montiel olea (which works with vec(At') opposed to vec(At)
n3 = n+1;
D = duplication(n+1); 
Dk = (D'*D)\D';  
L = elimination(n3);
Knn = commutation(n3,n3);     
D2 = duplication(2*(n+1)); 
Dk2 = (D2'*D2)\D2';   
RForm.n = n;   RForm.p = p;  RForm.Tmp = T-p;  RForm.H = H; 
if nargin<10  %% Select t_b if not provided
   a = 1; WaldStatsPM = zeros(length(idx),1);
   for t = idx 
       %% VAR augmented with IV
       wt = W(:,a);
       wX3 = wt.*X3;
       At3 =  (wX3'*X3)\(wX3'*y3); % VAR-augmented
       Uhat3 = y3 - X3*At3;
       Sighat3small = 1/H*(Uhat3.*wt)'*(Uhat3);
       Sigg3 = 1/H*(Uhat3.*wt)'*(Uhat3.*wt);
       COVhatsigmasmall = 2*Dk*kron(Sigg3,Sighat3small)*Dk'; 
       Pchol = chol(Sighat3small,'lower');
       Hder = L'/(  L * (speye(n3^2)+Knn) * kron(Pchol,speye(n3)) * L'  );
       VarChol = Hder*COVhatsigmasmall*Hder';
       WaldStatsPM(a) = (((H^.5)*Pchol(1+idx_shock))^2)/VarChol(1+idx_shock,1+idx_shock);
       a = a + 1; 
   end 
   [~,t_stand_idx]=max(WaldStatsPM);
   t_standard = idx(t_stand_idx); 
end
output.t_stand = t_standard; 
output.t_stand_idx = find(idx==t_standard);    




wt_b =  H*GaussK((1:T-p)',t_standard,H)./sum(GaussK((1:T-p)',t_standard,H),1);    

%% VAR augmented with IV
wtbX3 = wt_b.*X3;
Abt3 =  (wtbX3'*X3)\(wtbX3'*y3); % VAR-augmented
Uhatb3 = y3 - X3*Abt3;
indicesAraw = reshape(1:n3*k3,k3,n3)'; 
indicesA_3 = indicesAraw(:,ic+1:end);
reidxPMvecA = vec(indicesA_3); 
IDXsig = tril(ones(2*n3)); nsig = sum(sum(IDXsig));
IDXsig(IDXsig==1)=1:nsig;
El1 = IDXsig(1:n3,1:n3); El1 = vec(El1(El1~=0));
El2 = IDXsig(n3+1:2*n3,n3+1:2*n3); El2 = vec(El2(El2~=0));
Isig1 = eye(nsig); Isig1 = Isig1(El1,:);
Isig2 = eye(nsig); Isig2 = Isig2(El2,:);



%% Comopute IRFs at each t in "idx"    
output.VARpm = zeros(horizons+1,n,length(idx));   
output.VARpm_lb_msw = zeros(horizons+1,n,length(idx));
output.VARpm_ub_msw = zeros(horizons+1,n,length(idx));  
output.VARpm_lb_dm = zeros(horizons+1,n,length(idx));
output.VARpm_ub_dm = zeros(horizons+1,n,length(idx));

a = 1;
for t = idx    
    wt = W(:,a);    
    %% VAR augmented with IV 
    wX3 = wt.*X3;
    At3 =  (wX3'*X3)\(wX3'*y3); % VAR-augmented 
    Uhat3 = y3 - X3*At3; 
    Sighat3small = 1/H*(Uhat3.*wt)'*(Uhat3);   
    Sigg3 = 1/H*(Uhat3.*wt)'*(Uhat3.*wt); 
    COVhatsigmasmall = 2*Dk*kron(Sigg3,Sighat3small)*Dk';   
    iPix = inv(1./H*wX3'*X3); Piww =1./H*wX3'*wX3; 
    COVhat3 = kron(Sighat3small,iPix*Piww*iPix);  %*Piww*iPix
    Pchol = chol(Sighat3small,'lower');      
    Hder = L'/(  L * (speye(n3^2)+Knn) * kron(Pchol,speye(n3)) * L'  ); 
    VarChol = Hder*COVhatsigmasmall*Hder';      
    output.WaldStatPM(a) = (((H^.5)*Pchol(1+idx_shock))^2)/VarChol(1+idx_shock,1+idx_shock); 
    % Compute joint distribution of errors 
    wtall2 = [wt.*ones(1,n+1), wt_b.*ones(1,n+1)] ;     
    UhatAll2 = [Uhat3, Uhatb3]; 
    UhatAll2w =  UhatAll2.*wtall2 ;
    Sighat3 = 1/H*((wtall2.*UhatAll2)'*UhatAll2);      
    Sigg3 = 1/H*(UhatAll2w)'*(UhatAll2w);   
    COVhatsigma = 2*Dk2*kron(Sigg3,Sighat3)*Dk2';     
    Pchol2a = chol(Sighat3(1:n3,1:n3),'lower');
    Pchol2b = chol(Sighat3(n3+1:2*n3,n3+1:2*n3),'lower'); 
    HderA = (L'/(  L * (speye((n3)^2)+Knn) * kron(Pchol2a,speye(n3)) * L'  ))*Isig1; 
    HderB = (L'/(  L * (speye((n3)^2)+Knn) * kron(Pchol2b,speye(n3)) * L'  ))*Isig2; 
    Hder = [HderA;HderB];
    vPchol = [vec(Pchol2a); vec(Pchol2b)];
    VarChols = Hder*COVhatsigma*Hder';    
    idxchol = [1:n3,(n3^2)+idx_shock+1]; % plus one since you're looking at the IV augmented model
    RForm.WHat = blkdiag(COVhat3(reidxPMvecA,reidxPMvecA),... 
        VarChols(idxchol,idxchol));  
    RForm.AL = At3(ic+1:end,:)'; 
    RForm.Gamma  = vPchol(idxchol);  
    % Compute IRFs based on Asymptotic Distribution:
    % sqrt(T)*vec(Th_i-Th) ->
    % N(0, C_i AValpha C_i' + Cb_i AVbeta Cb_i' +
    %           C_i*ACOValb*Cb_i' + Cb_i*ACOValb*C_i') 
    AA = [RForm.AL;[speye(n3*(p-1)),sparse(n3*(p-1),n3)]];
    JJ = [speye(n3),sparse(n3,n3*(p-1))]; 
    Apoweri = AA;
    Phis = zeros(n3,n3,horizons+1); 
    Phis(:,:,1) = eye(n3);
    for i = 1:horizons
        Phis(:,:,i+1) = JJ*Apoweri*JJ';
        Apoweri = Apoweri * AA;
    end 
    I1 = eye(n3+1); I1(end,:)=[]; I2 = eye(n3+1); I2(1:end-1,:)=[]; 
    B_it_lam2 = (I1*RForm.Gamma)./(I2*RForm.Gamma);  
    dBdGam_2 = I1./(I2*RForm.Gamma) - [zeros(n3),(I1*RForm.Gamma)./((I2*RForm.Gamma)^2)];  
    Theta_2 = zeros(n3,horizons+1);
    Theta_var_dm_2 = zeros(n3,horizons+1); 
    Theta_lb_AR_2 = zeros(n3,horizons+1);
    Theta_ub_AR_2 = zeros(n3,horizons+1);
    for i = 0:horizons
        Theta_2(:, i+1) = Phis(:,:,i+1)*B_it_lam2;
        G_i = zeros(n3^2,n3^2*p);
        for mm = 0:i-1
            G_i = G_i + kron(JJ*(AA')^(i-1-mm),Phis(:,:,mm+1)) ;
        end
        C_i = ( kron(B_it_lam2',speye(n3))*G_i );
        Cbar_i =  Phis(:,:,i+1)*dBdGam_2 ; 
        GradAG = [C_i, Cbar_i]; 
        COVthetat = GradAG*RForm.WHat*GradAG'; 
        var_diag = diag( COVthetat );
        var_diag(var_diag<0)=0;
        Theta_var_dm_2(:,i+1) = var_diag;
        
        % AR CI
        B_it2 =  I1*RForm.Gamma  ;   
        C_i = ( kron(B_it2',speye(n3))*G_i );
        Cbar_i2 =  Phis(:,:,i+1)*I1 ; 
        GradAG_AR2 = [[C_i, Cbar_i2  ];[zeros(1,size(GradAG,2)-1),1]]; 
        COVthetat_AR2 = GradAG_AR2*RForm.WHat*GradAG_AR2'; 
        %full(RForm.WHat(end-4:end,end-4:end))
        e = eye(n3);  ep1 = eye(n3+1);
        critval = norminv(1-alpha/2,0,1)^2;
        CIs2 = zeros(2,n3);
        for ii = 1:n3
              XYf =  [e(ii,:)*Phis(:,:,i+1)*B_it2; I2*RForm.Gamma];
              What = ep1([ii,n3+1],:)*COVthetat_AR2*ep1([ii,n3+1],:)'/H;  
              f0 = XYf(1)^2-critval*What(1,1);
              f1 = XYf(1)*XYf(2) - critval*What(1,2);
              f2 = XYf(2)^2 - critval*What(2,2);
              D = f1^2 - f0*f2;
              if D<0
                  D = 0;
              end
              r1 = (f1-D^(1/2))/f2;
              r2 = (f1+D^(1/2))/f2;
              if and(D>=0,f2>=0)
                  CIs2(:,ii) = [r1,r2];
              elseif and(D<0,f2<0)
                  CIs2(:,ii) = [-inf,inf];
              elseif and(D>=0,f2<0)
                 CIs2(:,ii) =  [-inf,inf];
              else
                  disp('stop')
              end
        end
        Theta_lb_AR_2(:,i+1) = CIs2(1,:)';
        Theta_ub_AR_2(:,i+1) = CIs2(2,:)'; 
    end
    SN_q = norminv(1-alpha/2,0,1);
    output.VARpm(:,:,a) = Theta_2(2:end,:)';
    IRF_sd = sqrt( Theta_var_dm_2 ./H ); 
    output.VARpm_ub_msw(:,:,a) = Theta_ub_AR_2(2:end,:)';
    output.VARpm_lb_msw(:,:,a) = Theta_lb_AR_2(2:end,:)';
    output.VARpm_ub_dm(:,:,a) = (Theta_2(2:end,:) + IRF_sd(2:end,:).*SN_q)';
    output.VARpm_lb_dm(:,:,a) = (Theta_2(2:end,:) - IRF_sd(2:end,:).*SN_q)';
    
    
    
    a = a + 1;
end
 

end

 