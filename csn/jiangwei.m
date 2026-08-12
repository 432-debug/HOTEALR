% 加载参考数据
data1 = load('E:\cd45NOR_ICA\jiangwei\GIN_network_information.mat');
data2 = load('E:\cd45NOR_ICA\mouse_scRNAseq_data.mat');

% 定义输入路径前缀和输出目录
path_in_prefix = 'E:\cd45NOR_ICA\mouse_scRNAseq_net\26W_ICA_sub_network5\mouse_scRNAseq_network_';
out_dir = 'E:\cd45NOR_ICA\mouse_scRNAseq_net_jiangwei\26W_ICA_sub_network5\';
path_out_prefix = [out_dir, 'mouse_scRNAseq_jiangwei_'];

% 如果输出文件夹不存在，则自动创建
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% 处理标签，获取需要保留的索引
labels = data2.Final_genes; 
matrixData = data1.reliable_interactions;
col1 = matrixData(:, 1); 
unique_labels = unique(col1, 'stable'); 

% 获取相同标签在原 Final_genes 中的索引 (idx_first)
[same_labels, ~, idx_first] = intersect(unique_labels, labels, 'stable');
idx_first_matrix = idx_first(:); % 转为列向量
idx_first_sorted = sort(idx_first_matrix, 'ascend'); % 升序排序，保持原基因相对顺序

% 开始循环降维
for fnum = 1:20
    % 1. 正确拼接完整的文件路径 (使用中括号拼接字符串)
    file_in = [path_in_prefix, num2str(fnum), '.mat'];
    file_out = [path_out_prefix, num2str(fnum), '.mat'];
    
    % 2. 加载原始邻接矩阵
    data3 = load(file_in);
    test_adjacency = data3.adjacency;
    
    % 3. 执行降维 (切片)
    % 注意：这里直接把降维后的矩阵命名为 adjacency，这样 save 时变量名才能对得上
    adjacency = test_adjacency(idx_first_sorted, idx_first_sorted);
    
    % 4. 保存为新的 .mat 文件
    save(file_out, 'adjacency');
    
    % 打印进度，方便观察
    fprintf('已完成: %s (降维后维度: %d x %d)\n', file_out, length(idx_first_sorted), length(idx_first_sorted));
end

disp('全部 20 个矩阵降维完成！');