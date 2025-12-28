function sum = submatrix_handle(A_J,b_j,x)
%SUBMATRIX_HANDLE 输入矩阵A_j,b_j当前解x，计算每次kaczmarz迭代的加和
%   此处显示详细说明

m = size(A_J,1);
n = size(x,1);
norm_rows = vecnorm(A_J,2,2);
%默认权重为 1/m
sum = zeros(n,1);
w = 1/m;
for i = 1:m
    sum = sum +  w *(A_J(i,:)*x -b_j(i) )/norm_rows(i)^2*A_J(i,:)';
end


end