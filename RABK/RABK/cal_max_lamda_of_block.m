function lamda = cal_max_lamda_of_block(A,J)
%MAX_LAMDA_OF_BLOCK 
%求所有块矩阵的最大lamda
% 输入 矩阵A ，分块索引J，

num_block = size(J,1);

lst = zeros(num_block,1);

 %%求AJ' * diag(1/||  ai  ||^2) * AJ 的特征值,并对所有J取最大

for i = 1 : num_block
    rows_used = J{i} ; 
    A_J = A(rows_used,:);
    rows_norm = vecnorm(A_J,2,2);
    diag_vec = 1./(rows_norm.^2) ;
    D = diag(diag_vec);
    W = A_J' * D * A_J ;
    char_value = eig(W);
    lamda_MAX = max(char_value); % 每个块的最大lamda
    lst(i) = lamda_MAX;
end


lamda = max(lst);


end
