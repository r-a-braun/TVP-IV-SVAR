function [IP] = compute_special_IP(X,RIW)
%COMPUTE_IP Summary of this function goes here
%   Detailed explanation goes here  
T = length(X);
IP = ones(T,1);
Xgrow = X(2:end,:)./X(1:end-1,:);
RIWlag = RIW(1:end-1,:); 
for t = 2:T 
    IP(t) = IP(t-1)*( nansum(RIWlag(t-1,:).*Xgrow(t-1,:)) ./ nansum(RIWlag(t-1,:) ) ); 
end

