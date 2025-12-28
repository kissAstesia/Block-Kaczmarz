function J = row_partition(A,t)
%行分块函数，输入A，块的尺寸t，输出所有块的行索引
m = size(A,1);
num_block = ceil(m/t);
J = cell(num_block,1);

    for i = 1:num_block
        start_idx = 1 + (i-1)*t;
        end_idx = min(i*t,m);
        J{i} = start_idx:end_idx;
    
    end
end