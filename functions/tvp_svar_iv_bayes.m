function output = tvp_svar_iv_bayes( Yraw, ziv, p ,  X_exo, idx, horizons, idx_stationary)
% Function to compute Bayesian estimator
 
Yraw(:,idx_stationary==0) = [ NaN(1, sum(idx_stationary==0)); diff(Yraw(:,idx_stationary==0)) ];
if sum(idx_stationary==0)>0
    Yraw = Yraw(2:end,:); ziv = ziv(2:end,:); X_exo = X_exo(2:end,:); idx = idx - 1;
end

[TpP , n] = size(Yraw);
T = TpP-p;
zp = ziv(p+1:end);
Lags = lagmatrix(Yraw,1:p); 
X_exo = X_exo(p+1:end,:);
X = [X_exo, zp, ones(T,1), Lags(p+1:end,:)]';
y = Yraw(p+1:end,:)';   
Ahat = (y*X')/(X*X');
dummies_hat = Ahat(:,1:size(X_exo,2))';
Ahat_zX = Ahat(:, size(X_exo,2)+1:end)';   
Uhat = y - Ahat*X;  
Sighat = cov(Uhat');

%% Preliminaries 
shortY = Yraw(p+1:end,:); 
Y = reshape(shortY',T*n,1); 
X = lagmatrix(Yraw,1:p); 
X = X(p+1:end,:);  
X_all = [zp, ones(T,1), X]; 
n_param = size(X_all,2)*n; 
tempid = reshape(1:T*n_param,n_param,T)';
idi1 = kron((1:n*T)',ones(size(X_all,2),1));
idj1 = reshape(tempid(:,1:n_param)',T*n_param,1); 
bigX = sparse( idi1 , idj1 , reshape(kron(X_all,ones(n,1))',T*n_param,1) );  
k  = size(X_all,2)*n    ; % dimension of states 
Htheta = speye(T*k) - sparse(k+1:T*k,1:(T-1)*k,ones((T-1)*k,1),T*k,T*k);
 
 
%% prior 
atheta = zeros(k,1); Vtheta = 10*ones(k,1);
nutheta0 = 5*ones(k,1); 
Stheta0 = .01^2*ones(k,1).*(nutheta0-1);      
S0sigma = eye(n); nu0sigma = n+3; 
nsims = 10000; burnin = 1000;

%% initialize 
Sigtheta = .01*ones(k,1);
Sig = Sighat;
theta0 = vec(Ahat_zX);
n_exo = length(dummies_hat(:)) ;
%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting TVP-VARX.... ');
start_time = clock;
store_irf = zeros(horizons+1,n,length(idx),nsims);
        JJ = [speye(n),sparse(n,n*(p-1))];  

for isim = 1:nsims+burnin    
   
     %% sample TVP parameters  
    ytil = Y - vec((X_exo*dummies_hat)');
    invS = sparse(1:T*k,1:T*k,repmat(1./Sigtheta',1,T),T*k,T*k);
    invSig = kron(speye(T),inv(Sig)); 
    XinvSig = bigX'*invSig;
    HinvSH = Htheta'*invS*Htheta;
    alptheta = Htheta\[theta0;sparse((T-1)*k,1)];   
    Ktheta = HinvSH + XinvSig*bigX;
    dtheta = XinvSig*ytil + HinvSH*alptheta;
    thetahat = Ktheta\dtheta;
    theta = thetahat + chol(Ktheta,'lower')'\randn(T*k,1); 
    Theta = reshape(theta,k,T)';
    uhat = ytil - bigX*theta; 
    Uhat = reshape(uhat,n,T)';
    SSE = Uhat'*Uhat;
    Sig = iwishrnd(SSE + S0sigma,nu0sigma+T);   
    
    %% Sample the dummies (if there are)
    Ytil_dum = reshape(Y - bigX*theta,n,T)'; 
    Z  = kron(eye(n),X_exo);
    VARIANCE = kron(inv(Sig),speye(T));
    V_post_dummies = inv( Z'*VARIANCE*Z);  
    a_post_dummies = V_post_dummies*( Z'*VARIANCE*Ytil_dum(:) ); 
    vec_dummies =  a_post_dummies + chol(V_post_dummies)'*randn(n_exo,1); 
    dummies_hat = reshape(vec_dummies,size(dummies_hat,1),size(dummies_hat,2));

    %% sample theta0
    Ktheta0 = sparse(1:k,1:k,1./Sigtheta + 1./Vtheta);
    theta0hat = Ktheta0\(atheta./Vtheta + theta(1:k)./Sigtheta);
    theta0 = theta0hat + chol(Ktheta0,'lower')'\randn(k,1);
          
   
    %% sample Sigtheta
    e = reshape(theta-[theta0;theta(1:(T-1)*k)],k,T);
    Sigtheta = 1./gamrnd(nutheta0+T/2, 1./(Stheta0 + sum(e.^2,2)/2));       


    if isim>burnin 

        %% Compute IRFs at different time poinds:
        IRFs_idx = zeros(horizons+1,n,length(idx));
        Theta_idx = Theta(idx,:);
        Theta_1 = reshape(Theta_idx(1,:)',k/n,n);
        for j = 1:size(Theta_idx,1) 
            Theta_j = reshape(Theta_idx(j,:)',k/n,n);
            irf_rpo_h0t = (Theta_j(1,:)./Theta_1(1,1))';
            AA = [Theta_j(3:end,:)'; [speye(n*(p-1)),sparse(n*(p-1),n)]]; 
            Phis = zeros(n,horizons+1);
            Phis(:,1) = irf_rpo_h0t;
            Apoweri = AA;
            for ii = 1:horizons
                Phis(:,ii+1) = JJ*Apoweri*JJ'*irf_rpo_h0t;
                Apoweri = Apoweri * AA;
            end
            IRFs_j = Phis' ;
            IRFs_j(:,idx_stationary==0)=cumsum(IRFs_j(:,idx_stationary==0));
            IRFs_idx(:,:,j) = IRFs_j;
        end

        i = isim-burnin;
        store_irf(:,:,:,i) =  IRFs_idx;      
    end
    
    if ( mod( isim, 10000 ) ==0 )
        disp(  [ num2str( isim ) ' loops... ' ] )
    end 
    if ( mod( isim, 10000 ) ==0 )
        disp(  [ num2str( isim ) ' loops... ' ] )
    end 
    
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );


output.store_irf = store_irf; 


 

end
 