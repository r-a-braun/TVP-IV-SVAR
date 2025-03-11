function [lambda_Wald, lambda_F, pval_Wald, pval_F , pval_F_2, B, variance_Cbeta] = Granger_Test(y, p, X_exo, h_sel )
%Input: y: a (T+p)xK Matrix of Observations
%       p: lag order
%       intercept: a dummy which indicates, whether the model includes an
%       intercept
%       xz_vector: a vector, which indicates, which variable belongs to xt
%       and zt (xt:=1,zt:=0,else=2)


[B, Sigma_U, ~, ~ , Z] = VARls(y,p,X_exo);


vec_B = B(:); 
[T,K] = size(y); 

C = h_sel;
N = size(h_sel,1);
variance_Cbeta = C*kron(inv(Z*Z'),Sigma_U)*C'; 
lambda_Wald=(C*vec_B)'*inv(C*kron(inv(Z*Z'),Sigma_U)*C')*C*vec_B; % Test statistic for Wald test
lambda_F = lambda_Wald/N; % Corresponding F-Statistic


pval_Wald = 1-chi2cdf(lambda_Wald, N); % computing p value
pval_F = 1-fcdf(lambda_F, N, K*T-K^2*p-K); % computing p value
pval_F_2 = 1-fcdf(lambda_F, N,  T-K*p-1); % computing p value


end




function [Bhat, Sigmahat, Uhat, VarBhat , Z] = VARls(y, p, X_exo)
%VAR Computes LS estimates of VAR(p) with intercept
% inputs:   -y: the data of size (T+pxK)
%           -p: VAR lag order,
%           -int: intercept dummy
% outputs:  -Bhat: OLS estimates AR parameters
%           -Sigmahat: OLS estimate COV matrix
%           -Uhat: residuals
%           -Tstatistic:  t stats of parameters
[Traw, K] = size(y);
T = Traw - p;
Y = y(p+1:end,:)';  
rhs = lagmatrix(y,1:p);
Z = rhs(p+1:end,:)'; 
Z = [ X_exo(p+1:end,:)'; Z ]; 
Bhat = Y*Z'*inv(Z*Z'); %LS Estimator
Uhat = (Y-Bhat*Z);
Sigmahat = 1/(T-K*p-1)*Uhat*Uhat'; % Residual Cov-Matrix
VarBhat = kron(inv(Z*Z'),Sigmahat); % Cov Matrix of Estimated Parameters 
end


