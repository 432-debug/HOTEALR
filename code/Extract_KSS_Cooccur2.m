function KeySubgraphSet = Extract_KSS_Cooccur2(All_F1_Masks, test_adjacency, D_score)
    % 基于共现矩阵快速提取 KSS 节点集合，采用连接权重代理分数和子图分解。
    
    if isempty(All_F1_Masks)
        KeySubgraphSet = {};
        return;
    end
    
    [N_F1, D] = size(All_F1_Masks);
    
    % --- 0. 预计算全局平均绝对连接权重 ---
    % 仅考虑 test_adjacency 中非零的权重
    non_zero_weights = nonzeros(test_adjacency);
    absw = mean(abs(non_zero_weights)); % 全局平均绝对连接权重
    if absw == 0
        absw = 1e-6; % 避免除以零或阈值为零
    end
    
    % --- 1. 参数设置 & 共现矩阵 (不变) ---
    tau_cooccur = ceil(N_F1 * 0.4); 
    
    % --- 2. 快速计算共现矩阵 ---
    CoMatrix = All_F1_Masks' * All_F1_Masks;
    
    % --- 3. 提取高频共现边 (KSS 种子) ---
    [row, col] = find(CoMatrix >= tau_cooccur & triu(ones(D, D), 1));
    KSS_Seeds = unique([row, col], 'rows');
    
    if isempty(KSS_Seeds)
        KeySubgraphSet = {};
        return;
    end
    
    % --- 4. 子图分解、代理分数计算与筛选 ---
    KeySubgraphSet = {};
    processed_nodes = false(1, D); 
    
    % 代理分数阈值：初始化为保守值，后续可以动态调整
    proxy_threshold = 0.01; 

    for i = 1:size(KSS_Seeds, 1)
        node1 = KSS_Seeds(i, 1);
        node2 = KSS_Seeds(i, 2);
        
        % 如果种子节点已处理，跳过
        if processed_nodes(node1) && processed_nodes(node2)
            continue;
        end
        
        Ik_nodes = unique([node1, node2]);
        
        % --- 4.1 局部贪婪扩展 (保持原逻辑，基于 D_score 引导) ---
        % 这步基于 D_score 扩展是为了确保 KSS 拥有高重要性节点
        max_expand_steps = randi([1, 3]); 
        for step = 1:max_expand_steps
            current_D_scores = D_score(Ik_nodes);
            mean_D = mean(current_D_scores);
            
            neighbors = find(any(test_adjacency(Ik_nodes, :) ~= 0, 1));
            expand_candidates = setdiff(neighbors, Ik_nodes);
            
            if isempty(expand_candidates)
                break;
            end
            
            D_candidate = D_score(expand_candidates);
            [~, max_idx] = max(D_candidate);
            best_candidate = expand_candidates(max_idx);
            
            if D_score(best_candidate) > mean_D * 0.9 
                Ik_nodes = [Ik_nodes, best_candidate];
            else
                break;
            end
        end
%% “best_candidate = expand_candidates(max_idx);”改成 选取与现有节点Ik_nodes（expand_candidates）相连的节点中权重绝对值最大的       

        % --- 4.2 KSS 分解：基于平均连接权重 absw ---
        % 目的：将 Ik_nodes 分解为多个内部连接更紧密的小子图
        
        % 提取 Ik_nodes 内部的连接子图 (邻接矩阵)
        Ik_adj = test_adjacency(Ik_nodes, Ik_nodes);
        
        % 寻找弱连接：权重绝对值小于全局平均 absw 的边
        weak_links = (abs(Ik_adj) < absw);
        
        % 移除弱连接：构建强连接子图 Ik_strong_adj
        Ik_strong_adj = Ik_adj;
        Ik_strong_adj(weak_links) = 0;
        
        % 连通分量分析 (将强连接子图分解为子模块)
        % 假设您有一个 Graph/DFS/BFS 相关的函数来找到连通分量 (CC)
        % 如果没有现成的函数，可以使用简化方法或 MATLAB 的 graph/conncomp
        
        % 简单近似：对 Ik_strong_adj 进行连通分量分析
        CC = Weak_Link_Component_Analysis(Ik_strong_adj); % 假设存在此函数
        
        % 4.3 评估分解后的子图模块
        original_nodes_map = Ik_nodes;
        
        for k = 1:length(CC)
            cc_indices = CC{k}; % 在 Ik_nodes 内部的索引
            KSS_module_nodes = original_nodes_map(cc_indices); % 实际节点 ID
            
            if length(KSS_module_nodes) < 2
                continue; % 忽略单个节点的模块
            end
            
            % --- 4.4 KSS 代理分数修正：使用连接权重 ---
            % 提取 KSS 模块内部的连接权重矩阵
            Module_adj = test_adjacency(KSS_module_nodes, KSS_module_nodes);
            Module_weights = nonzeros(triu(Module_adj, 1)); % 上三角非零权重
            
            if isempty(Module_weights)
                proxy_score = 0; % 无内部连接，分数设为 0
            else
                % KSS Score = |内部连接权重|均值 * |内部连接权重|方差
                mean_W = mean(abs(Module_weights));
                std_W = std(abs(Module_weights), 0);
                proxy_score = mean_W * std_W;
            end
            
            % --- 4.5 最终筛选与存储 ---
            % 这里可以动态调整 proxy_threshold，但为了简化，我们先用一个保守值。
            if proxy_score >= proxy_threshold
                KSS_entry.Nodes = KSS_module_nodes;
                KSS_entry.Score = proxy_score;
                KeySubgraphSet{end+1} = KSS_entry;
                processed_nodes(KSS_module_nodes) = true; % 标记节点
            end
        end
    end
end
% 假设 Ik_strong_adj 是一个 NxN 的对称矩阵
function CC = Weak_Link_Component_Analysis(Ik_strong_adj)
    % 转换为布尔连通矩阵
    adj_bool = (Ik_strong_adj ~= 0);
    
    % 1. 创建图对象 (需要 Graph Theory Toolbox 或较新版本)
    G = graph(adj_bool); 
    
    % 2. 找到连通分量
    [bin, binsize] = conncomp(G); 
    
    CC = cell(length(binsize), 1);
    for k = 1:length(binsize)
        % 找到属于第 k 个分量的节点（在 Ik_nodes 中的局部索引）
        CC{k} = find(bin == k); 
    end
end