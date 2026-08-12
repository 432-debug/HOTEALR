function current_subspace = dynamic_subspace_pareto_neighbors(PopMask, FrontNo, edge_net)

    % 1. 取所有非支配前沿的个体（FrontNo == 1）
    F1_mask = PopMask(FrontNo == 1, :);
    
    % 若没有 Front 1 个体（极端情况），直接用全部种群
    if isempty(F1_mask)
        F1_mask = PopMask;
    end

    % 2. 找到 Front 1 个体用到的所有节点索引
    selected_nodes = find(any(F1_mask == 1, 1));   % 列维度上 OR

    % 3. 对每个节点取邻接节点
    neighbor_nodes = find(any(edge_net(selected_nodes, :) ~= 0, 1));

    % 4. 合并 + 去重
    current_subspace = unique([selected_nodes, neighbor_nodes]);

end
