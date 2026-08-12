function Offspring = Operator_KSS_Guided2(Parent, KSS, D_score, test_adjacency)
    % KSS 引导的结构变异算子
    % Parent 是 N x D 的矩阵
    [N_parent, D] = size(Parent);
    Offspring = zeros(N_parent, D); % 初始化 Offspring 矩阵
    
    % --- 1. 参数设置：激进增大子图规模 ---
    P_KSS_mutation = 0.4; % 40% 概率触发 KSS 结构变异
    P_KSS_add = 0.75;    % 75% 概率进行 KSS 添加 (解决规模小的问题)
    
    % 检查 KSS 是否为空
    if isempty(KSS)
        KSS_available = false;
    else
        KSS_available = true;
        % 预计算 KSS 概率
        KSS_scores = cellfun(@(x) x.Score, KSS);
        KSS_probs = KSS_scores / sum(KSS_scores);
    end
    % -------------------------------------------------------------
    % ★★★ 核心修正：循环处理每一个父代个体 ★★★
    % -------------------------------------------------------------
    for i = 1:N_parent
        current_parent = Parent(i, :);
        current_offspring = current_parent; % 从父代初始化后代
        
        % 检查是否触发 KSS 结构变异
        if KSS_available && rand() < P_KSS_mutation
            
            % 2. KSS 结构变异触发
            % 2.1 基于 KSS Score 加权选择一个 KSS 
            k_idx = randsample(length(KSS), 1, true, KSS_probs);
            Ik_nodes = KSS{k_idx}.Nodes;
            
            % 2.2 判断与当前父代的关系 (现在 current_parent 是 1x D 向量，all() 返回标量)
            % 检查 KSS 的所有节点是否都在当前父代中被激活
            is_contained = all(current_parent(Ik_nodes)); 
            if rand() < P_KSS_add % KSS 添加策略 (高概率)
                if ~is_contained 
                    % 添加操作：将 KSS 中的节点加入 Offspring
                    current_offspring(Ik_nodes) = true;
                else
                    % 如果已包含，执行常规变异以探索局部邻域
                    % 注意：Muation_T 必须能处理 1x D 向量
                    cand_off = Muation_T(current_offspring, D_score, test_adjacency, false);
                    if isempty(cand_off)
                        % 极端兜底
                        current_offspring = current_parent;
                    elseif size(cand_off,1) == 1
                        current_offspring = cand_off;
                    else
                        % 可以随机 / 最小度 / 最大连通度
                        idx = randi(size(cand_off,1));
                        current_offspring = cand_off(idx,:);
                    end
                end         
            else % KSS 移除策略 (低概率)
                % 移除条件：KSS 必须被包含 && 移除后解不为空 (sum(current_parent) > length(Ik_nodes))
                % 现在两个操作数都是标量，可以使用 &&
                if is_contained && sum(current_parent) > length(Ik_nodes) 
                    % 移除操作：将 KSS 中的节点移除 Offspring
                    current_offspring(Ik_nodes) = false;
                else
                    % 如果无法移除，执行常规变异
                    cand_off = Muation_T(current_offspring, D_score, test_adjacency, false);
                    if isempty(cand_off)
                        % 极端兜底
                        current_offspring = current_parent;
                    elseif size(cand_off,1) == 1
                        current_offspring = cand_off;
                    else
                        % 可以随机 / 最小度 / 最大连通度
                        idx = randi(size(cand_off,1));
                        current_offspring = cand_off(idx,:);
                    end  
                end
            end
        else
            % 触发常规变异 (低概率或 KSS 不可用时)
            % 传入当前个体的 Site 状态
            cand_off = Muation_T(current_offspring, D_score, test_adjacency, false);
            if isempty(cand_off)
                % 极端兜底
                current_offspring = current_parent;
            elseif size(cand_off,1) == 1
                current_offspring = cand_off;
            else
                % 可以随机 / 最小度 / 最大连通度
                idx = randi(size(cand_off,1));
                current_offspring = cand_off(idx,:);
            end
        end
        Offspring(i, :) = current_offspring;
    end
end