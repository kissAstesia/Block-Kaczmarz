function [x,iter] = RaBK(A,b,x_0,max_epoch,tol,t)
%RABK 随机块Kaczmarz方法，并行式，有外推步长
%应包含功能如下：
%1.输入系数矩阵A，向量b，初始解X_0，最大迭代次数Max_epoch,误差限tol，分块大小t
%默认使用partition sampling进行对子矩阵J 的取样
%输出迭代结束的解，以及迭代的次数


%为分析方便，若A的行数不可以整除块的尺寸t，则保留最后一个块，其权重w相应产生变化，其余保持不变


J = row_partition(A,t);  %得到块索引

num_block = size(J,1);



lamda_of_block = zeros(num_block,1);

for l = 1 : num_block
    [~ , lamda] = submatrix_handle(A(J{l},:),b(J{l}),x_0);
    lamda_of_block(l) = lamda;
end

max_lamda_of_block = max(lamda_of_block);          %提前得到lamda_max_block

iter = 0;                                          %迭代次数

%开始迭代，每次迭代会遍历一次所有的块，因此需要生成一个随机排列，大小为块的数量

    for i = 1: max_epoch
    
        perm = randperm(num_block);          %生成随机排列
    
        for j = 1 : num_block
            rows_used = J{perm(j)};
            A_J = A(rows_used,:);             %选取子矩阵
            b_j = b(rows_used);
            [sum,~] = submatrix_handle(A_J,b_j,x_0);


    
            alpha_k = 2*t/max_lamda_of_block;                        %%%在这里改变步长

            
    
            %进行迭代
            x_0 = x_0 - alpha_k * sum;
        end
    
    iter = iter + 1;
    
    if norm(A*x_0 - b,2) < tol
        fprintf('convergence at iteration %d',iter);
        break
    end
    
    end

x = x_0;
end