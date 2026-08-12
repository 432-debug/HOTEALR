function [subMasks, All_Masks, initial_subspace, W_fusion] = motif_vssps_enhanced(motif_label, chain_vec, test_adjacency, total_pop_size, D_score)
    % 目标：在初始化阶段激进地利用 D_score 和 chain_vec 生成高质量的Masks
    
    D = length(motif_label);
    clusters = unique(motif_label);
    K = length(clusters);
    % 强制转换为列向量，确保维度兼容
    D_score = double(D_score(:)); 
    chain_vec = double(chain_vec(:));
    % --- Step 1: 节点权重融合 (D_score 优先) ---
    % 归一化 D_score 和 chain_vec
    norm_D_score = (D_score - min(D_score)) / (max(D_score) - min(D_score) + eps);
    norm_chain_vec = (chain_vec - min(chain_vec)) / (max(chain_vec) - min(chain_vec) + eps);
    
    % W_fusion = 0.7 * D_score + 0.3 * chain_vec (赋予 D_score 更高的权重)
    % alpha = 0.7; 
    alpha = 1; 
    W_fusion = alpha * norm_D_score + (1 - alpha) * norm_chain_vec;
    
    % --- Step 2: 聚类大小 + 分配种群规模 (保持原逻辑) ---
    cluster_sizes = zeros(1, K);
    cluster_nodes = cell(1, K);
    for k = 1:K
        c_nodes = find(motif_label == clusters(k));
        cluster_nodes{k} = c_nodes;
        cluster_sizes(k) = length(c_nodes);
    end
    total_cluster_size = sum(cluster_sizes);
    sub_pop_sizes = round(total_pop_size * cluster_sizes / total_cluster_size);
    sub_pop_sizes(1) = sub_pop_sizes(1) + (total_pop_size - sum(sub_pop_sizes));
    
    % --- Step 3: 输出结构 ---
    All_Masks = false(total_pop_size, D);
    subMasks = cell(1, K);
    initial_subspace = cell(1, K); % 注意：这里返回的将是聚焦后的子空间
    global_idx = 1;
    min_len = 3;
    max_len = 5;
    
    % fprintf("=== 基于 motif 聚类 + 权重融合 (D_score 优先) 的激进链式初始化 ===\n");
    
    % --- Step 4: 每个聚类生成激进采样的子种群 ---
    for k = 1:K
        nodes_k = cluster_nodes{k};
        pop_k = sub_pop_sizes(k);
        mask_k = false(pop_k, D);
        
        W_fusion_k = W_fusion(nodes_k);
        
        % 找出该聚类内 W_fusion 最高的节点作为首选中心节点
        [~, max_idx_rel] = max(W_fusion_k);
        center_node_best = nodes_k(max_idx_rel); 

        fprintf("聚类 %d：节点 %d → 个体 %d，最高权重中心：%d\n", k, length(nodes_k), pop_k, center_node_best);
        
        % 策略 A: 激进采样 (前 1/3 个体围绕最高权重节点)
        num_aggressive = ceil(pop_k / 3); 
        
        for i = 1:pop_k
            % ------------------- 激进中心节点选择 -------------------
            if i <= num_aggressive
                % 激进利用：强制使用最高权重的节点作为中心
                center = center_node_best;
            else
                % 平衡探索：从聚类内随机选择中心节点，但概率偏向高权重
                P_selection = W_fusion_k / sum(W_fusion_k);
                center_idx_rel = randsample(length(nodes_k), 1, true, P_selection); % 基于概率采样
                center = nodes_k(center_idx_rel);
            end
            
            sel = center;
            s_len = randi([min_len, max_len]);
            current_node = center;
            
            % ------------------- 链式扩展（W_fusion 引导） -------------------
            while length(sel) < s_len
                neigh_D = find(test_adjacency(current_node, :) ~= 0);
                neigh_k = intersect(neigh_D, nodes_k);  % 只考虑聚类内节点
                neigh_k = neigh_k(~ismember(neigh_k, sel)); % 避免环路/重复
                
                if isempty(neigh_k)
                    break;
                end
                
                % 按 W_fusion 最大选邻居
                W_neigh = W_fusion(neigh_k);
                [~, max_idx] = max(W_neigh);
                next_node = neigh_k(max_idx);
                
                sel(end+1) = next_node;
                current_node = next_node;
            end
            
            % 写入 mask
            mask_k(i, sel) = true;
            All_Masks(global_idx, sel) = true;
            global_idx = global_idx + 1;
        end
        
        mask_k = unique(mask_k, 'rows', 'stable');
        subMasks{k} = mask_k;      

        % ------------------- Step 5: 初始子空间聚焦 (Top 30%) -------------------
        Top_X_percent = 0.3; % 仅保留该聚类内 W_fusion 最高的 Top 30% 节点
        [~, W_rank] = sort(W_fusion_k, 'descend');
        S_top_k = ceil(length(nodes_k) * Top_X_percent);
        
        % 返回 W_fusion 最高的 Top 30% 节点作为该子种群的初始子空间
        initial_subspace{k} = nodes_k(W_rank(1:S_top_k));
    end
    % fprintf("\n激进初始化完成：生成 %d 个子种群，总个体数 %d\n", K, sum(cellfun(@(x) size(x,1), subMasks)));
end
