function [H, CVs ] = choose_H_CF(Ytil, p,ic, Hgrid, nthin, hfc, var_choice, fraction_sample) 

[T , n] = size(Ytil);
if nargin < 6
    hfc = 1; 
    var_choice = 2:size(Ytil,2);
    fraction_sample = 1/2;
end
 
Lags = lagmatrix(Ytil,1:p);
Xtil = [ones(T-p,ic), Lags(p+1:end,:)];
ytil = Ytil(p+1:end,:);
k = n*p+ic;
GaussK = @(j,t,H) exp(-0.5*((j-t)/H).^2);
idxraw = (ceil((T-p)*(fraction_sample)):nthin:T)-p-hfc;


idx = intersect( find(ytil(:,1)~=0)-1,idxraw )' ;
warning('off')
CVs = NaN(1,length(Hgrid));
fbar = NaN(hfc,n); fbar(:,1) = 0; 
Cbar = eye(n*hfc);
Cbar(vec(isnan(fbar')==1),:)=[];
 
pvar = 13; 
weights_ivariance = zeros(1,n);
for i = 1:n
    yi = ytil(:,i);
    xi = [ones(length(yi),1) , lagmatrix(yi,1:pvar)];  yi = yi(pvar+1:end); xi = xi(pvar+1:end,:);
    beta = (xi'*xi)\(xi'*yi); ei  = yi-xi*beta;
    weights_ivariance(i) = 1./var(ei);
end
weights_ivariance = weights_ivariance./sum(weights_ivariance);

 
CVs2 = NaN(1,length(Hgrid));

for i = 1:length(Hgrid)
    %   clc;
    % disp(i/length(Hgrid))
    H = Hgrid(i);
    W =  H*GaussK((1:T-p)',idx,H)./sum(GaussK((1:T-p)',idx,H),1);
    MSEs = zeros(length(idx),1); 
    UhatIS = zeros(length(idx),n);
    
    for ii = 1:length(idx)
        ti = idx(ii);
        %% TVP model
        wt = W(:,ii);
        wt(ti+1:end)=0;
        wt = H.*(wt./sum(wt));
        wX = wt.*Xtil;
        At =  (wX'*Xtil)\(wX'*ytil); % VAR
        Uhat = ytil - Xtil*At;
        UhatIS(ii,:) = ytil(ti,:) - Xtil(ti,:)*At;
        Sighat_t = 1/H*(Uhat.*wt)'*(Uhat);  

        % Construct conditional forecast 
        fbar(:,1) = ytil(ti+1:ti+hfc,1);
        fbarvec = vec(fbar');
        fbarvec(isnan(fbarvec))=[];
        ytilfc = ytil( ti-p-1:ti,:)  ;  
        [bfc, Mfc] = get_fcmatrices(At',chol(Sighat_t,'lower'),hfc,ytilfc); 
        D = Cbar*Mfc';
        Dstar = pinv(D);
        muy = bfc + Mfc'*Dstar*(fbarvec - Cbar*bfc);
        Yhat = reshape(muy,n,hfc)';
        MSEs(ii) =  weights_ivariance(var_choice)*sum( ( ytil(ti+1:ti+hfc,var_choice) - Yhat(:,var_choice) ).^2 , 1)';   
    end 
    CVs(i) = sum(MSEs);  
end
[ ~ , idx_min] = min(CVs);
H = Hgrid(idx_min); 
end



