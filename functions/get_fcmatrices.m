function [b,M] = get_fcmatrices(Beta,P,h,Y)
% Constructs matrices as in Appendix A of scenario analysis paper
% Beta = reduced form matrix
% y_t = B*x + P* e_t     % REDUCED FORM 
% y_t'*A_0 = x_t'* A_p + e_t'   % structural form
% A_0' y_t = A_p' x_t + e_t      % structural form 2
% P^{-1}' --> A_0
% A_p = (P^{-1} B)' --> B*P^{-1}' --> B'*A_0  
[~,K]=size(Y);
p = (size(Beta,2)-1)/K;
[ Ac, Jc , nuc] = companion_internal( Beta , 1);
braw = zeros(h,K); % deterministic component
nuh = nuc;
YYt = vec(Y(end:-1:end-p+1,:)');
AA = Ac;
Mraw = zeros(K,K,h);
Mraw(:,:,1) = P';
for i =1:h
    braw(i,:) = Jc*(nuh + AA*YYt);
    nuh = nuc + Ac*nuh;
    if i<h
    Mraw(:,:,i+1) =(Jc*AA*Jc'*P)';
    end 
    AA = AA*Ac;
end
b = vec(braw');
M = zeros(h*K);
for j = 1:h
    M((j-1)*K+1:j*K,:) = [zeros(K,(j-1)*K),reshape(Mraw(:,:,1:h-j+1),K,K*(h-j+1))];
end




end

function [ A, J , nu] = companion_internal( Bhat , inc)
% Input: [c, A_1,...,A_p] of size (K, N*K+1)
% Output: sparse A, VAR(1) form
[K, Kpinc]=size(Bhat);
p = (Kpinc-inc)/K; 
nu = zeros(K*p,1);
if inc == 1
    nu(1:K) = Bhat(:,1);
end  
A = [Bhat(:,1+inc:end);[speye(K*(p-1)),sparse(K*(p-1),K)]]; 
J = [speye(K),sparse(K,K*(p-1))];   
end