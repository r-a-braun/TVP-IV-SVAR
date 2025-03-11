function SelectionM = VECH(K)
% Bertsche, Braun (2018): Identification of Structural Vector
% Autoregressions by Stochastic Volatility
%
% function to find the selection matrix for the vech operator
%
%% Input:
% - K (1 x 1): dimension
%
%% Output:
% - Selection(K(K+1)/2 x K^2): selection matrix

M = ones(K);
for i =2:K
    M(1:i-1,i) = 0; %make lower triangular
end
idx = find(M==1); li = length(idx); %find nonzero elements
SB = zeros(K^2,li);
for i = 1:li
    SB(idx(i),i) = 1;
end
SelectionM = pinv(SB);