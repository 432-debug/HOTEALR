function Mask = pathway_aware_vssps_node(data, pair_score, path2pairs)
    % 输入：
    %   data2：包含edge_net、gg等PEN数据
    %   pair_score：1125×1，基因对通路得分
    %   path2pairs：通路→基因对索引映射
    % 输出：
    %   Mask：N×1125矩阵，种群掩码（1表示选中该基因对）
    
    N = 300;  % 种群大小（可调整）
    n_pairs = size(data.subnetwork_adjacency,1);  % 1125
    Mask = zeros(N, n_pairs);
    
    % 1. 筛选高得分基因对（通路共现得分≥0.5）
    high_score_pairs = find(pair_score >= 0.5);
    if isempty(high_score_pairs)
        high_score_pairs = 1:n_pairs;  % 兜底：无高得分则全选
    end
    
    % 2. 条纹稀疏采样参数
    stripe_min_len = 3;   % 最小条纹长度（基因对数量）
    stripe_max_len = 8;   % 最大条纹长度
    high_ratio = 0.6;     % 60%条纹聚焦高得分基因对
    
    for i = 1:N
        % 随机生成条纹长度
        stripe_len = randi([stripe_min_len, stripe_max_len]);
        
        % 选择条纹区域：优先高得分基因对
        if rand < high_ratio
            % 从高得分基因对中随机选条纹
            if length(high_score_pairs) <= stripe_len
                stripe = high_score_pairs;
            else
                start_idx = randi(length(high_score_pairs) - stripe_len + 1);
                stripe = high_score_pairs(start_idx : start_idx + stripe_len - 1);
            end
        else
            % 随机选择基因对（保留探索性）
            start_idx = randi(n_pairs - stripe_len + 1);
            stripe = start_idx : start_idx + stripe_len - 1;
        end
        
        % 确保条纹在PEN网络中连通（利用edge_net邻接矩阵）
        if length(stripe) > 1
            stripe = ensure_connected(stripe, data.subnetwork_adjacency);
        end
        Mask(i, stripe) = 1;
    end
end

% 辅助函数：确保条纹内的基因对在PEN中连通（基于edge_net）
function stripe = ensure_connected(stripe, edge_net)
    % 检查当前条纹是否连通
    sub_net = edge_net(stripe, stripe);
    if isempty(sub_net) || sum(sub_net(:)) == 0
        % 若不连通，替换为与首个基因对相连的基因对
        first = stripe(1);
        neighbors = find(edge_net(first,:) ~= 0);
        if ~isempty(neighbors)
            stripe = [first, neighbors(randi(min(length(neighbors), length(stripe)-1)))];
        end
    end
    stripe = unique(stripe);  % 去重
end