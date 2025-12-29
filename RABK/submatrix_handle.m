function [sum , lamda] = submatrix_handle(A_J,b_j,x)
%SUBMATRIX_HANDLE 输入矩阵A_j,b_j当前解x，计算每次kaczmarz迭代的加和
%   输出max_lamda_of_block


m = size(A_J,1);
n = size(x,1);
norm_rows = vecnorm(A_J,2,2);

%默认权重为 1/m
sum = zeros(n,1);
w = 1/m;
for i = 1:m
    sum = sum +  w *(A_J(i,:)*x -b_j(i) )/norm_rows(i)^2*A_J(i,:)';
end

%%求AJ' * diag(1/||  ai  ||^2) * AJ 的特征值
D = zeros(m,1);
for j = 1 : m
    D(j) = 1/ norm_rows(j)^2;
end

W = A_J' * diag(D) * A_J ;
char_value = eig(W);
lamda = max(char_value);
    

end