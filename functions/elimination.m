function l = elimination(n)
% Constructs Elmination Matrix Lütkepohl p. 663 
l = zeros(n * (n + 1) / 2, n^2);
a = zeros(n, n);
for i = 1:n
    for j = 1:n
        a(i, j) = 1;
        l(:,(j - 1) * n + i) = vech(a);
        a(i, j) = 0;
    end
end  
end

function res = vech(A)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

% vech function that stacks the lower triangular portion of a n*n matrix as
% a n*(n+1)/2 vector. Its inverse is the math function

[~, n] = size(A);
res = zeros(n * (n + 1) / 2, 1);

tempbound = 0;
for i = 1:n
    res(tempbound + 1:tempbound + n - i + 1, 1) = A(i:n, i);
    tempbound = tempbound + (n - i) + 1;
end  

end