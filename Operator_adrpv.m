function [OffMask,poss_s_num] = Operator_adrpv(ParentMask,Site,rbm,allOne,other,D_score,test_adjacency,g,pv,UpdateRatio,current_subspace)
    % 入参说明：
    %   ParentMask - 全维度父代池（N_parent×D_full，原始父代掩码）
    %   current_subspace - 子空间维度索引（1×D）
    % 其他入参不变
    
    D = size(ParentMask, 2);
    S = length(current_subspace);  % 子空间维度数
    other = find(other);  % 将other逻辑数组转为索引数组（哪些列是目标列）
    other_sub = intersect(other, current_subspace);% other_sub：current_subspace中属于other的列（索引数组，列数=目标维度数）  
    other_nonsubspace = setdiff(other, other_sub);
    %实际上需要用交叉变异算子Operator改动的只有other部分(other_sub部分)

    % 父代拆分（全维度）
    %Parent1Mask = ParentMask(1:floor(end/2), :);
    %Parent2Mask = ParentMask(floor(end/2)+1:end, :);
    %Site1 = Site(:,1:floor(end/2));
    %Site2 = Site(:,floor(end/2)+1:end);
%% ========================
    [N, D] = size(ParentMask);
    half = floor(N / 2); % 确保父代成对

    % --- [关键修复 1] 强制父代分割严格对齐 (忽略多余的奇数行) ---
    Parent1Mask = ParentMask(1:half, :);
    Parent2Mask = ParentMask(half+1:2*half, :);

    % --- [关键修复 2] 强制 Site 对称 (忽略传入 Site 的后半段，强制 Site2=Site1) ---
    % 确保 Site 是行向量
    Site = reshape(Site, 1, []);
    % 截取前一半
    if length(Site) < half
        % 防御性代码：如果传入的 Site 不够长，补 false
        Site1 = [Site(1:end), false(1, half - length(Site))];
    else
        Site1 = Site(1:half);
    end
    % 强制镜像：父代2的行为必须与父代1完全一致
    Site2 = Site1;
%% ==========================
    % 重组 Site 以便 RBM 使用 (仅用于提取 RBM 输入数据)
    % 注意：这里的 Site_combined 长度为 2*n_half，与 Parent1Mask/Parent2Mask 对应
    Site = [Site1, Site2];

    OffMaskT1 = [];

    % =======================
    % ① RBM 分支
    % =======================
    if any(Site)
        if isempty(other_sub)
            other_sub = current_subspace;
        end
        % ------- 新增：维度检查 -------
        if isempty(rbm) || isempty(rbm.Weight)
            use_rbm = false;
        elseif size(ParentMask(Site,other_sub),2) ~= size(rbm.Weight,1)% 维度不一致 = 必然出错 → 回退
            fprintf('[WARNING] RBM skipped — dim mismatch: X(%d) vs W(%d)\n',...
                size(ParentMask(Site,other_sub),2), size(rbm.Weight,1));
            use_rbm = false;
        else
            use_rbm = true;
        end
        % ============== RBM 生成子代 ==============
        if use_rbm
            X = ParentMask(Site,other_sub);
            poss = rbm.reduce(X);   % safe now
            [OffTemp,parent1_indices] = Cross2(sum(Site),size(poss,2),poss);
            OffTemp = rbm.recover(OffTemp);

            OffMaskT1 = false(size(OffTemp,1), D);
            OffMaskT1(:,allOne) = true;
            OffMaskT1(:,other_sub) = OffTemp;

            ParentMask_Site = ParentMask(Site,:);
            OffMaskT1(:,other_nonsubspace) = ParentMask_Site(parent1_indices, other_nonsubspace);
        else
            % ========== fallback: 非 RBM pv 交叉 ==========
            % 相当于 treat 全部 Site 为非 RBM 个体
            Site1_f = Site(:,1:floor(end/2));
            Site2_f = Site(:,floor(end/2)+1:end);

            OffMaskT1 = Cross_pv( ...
                Parent1Mask(Site1_f,:), Parent2Mask(Site2_f,:), other, pv, current_subspace );
        end
    end

%    if any(Site)
%        if isempty(other_sub)
%            other_sub = current_subspace;
%        end 
%        disp(size(ParentMask(Site,other_sub)));
%        disp(size(rbm.Weight));
        
%        poss    = rbm.reduce(ParentMask(Site,other_sub));
%        [OffTemp,parent1_indices] = Cross1(sum(Site),size(poss,2),poss);
%        OffTemp = rbm.recover(OffTemp);

%        OffMaskT1 = false(size(OffTemp,1),size(Parent1Mask,2));
%        OffMaskT1(:,allOne) = true;
%        OffMaskT1(:,other_sub) = OffTemp;
%        ParentMask_Site = ParentMask(Site, :);
%        OffMaskT1(:,other_nonsubspace) = ParentMask_Site(parent1_indices, other_nonsubspace);
%    end

    %% 非RBM父代的Cross_pv逻辑（同理改造，确保非子空间继承Parent1Mask）
    OffMaskT2 = false(floor(sum(~Site)/2),size(ParentMask,2));
    if ~isempty(other_sub)

        [OffMaskT2] = Cross_pv(Parent1Mask(~Site1,:), Parent2Mask(~Site2,:), other, pv, current_subspace);
    end
    %OffMaskT2(:, allOne) = true;
    
    %% 合并子代+变异（Muation_T仅子空间变异）
    OffMask = [OffMaskT1; OffMaskT2];
    if isempty(OffMask)
        poss_s_num = 0;
        return;
    end
    [OffMask,poss_s_num] = Muation_T(OffMask,D_score,test_adjacency,Site);
end

function pop=creatpop(popnum,D,D_score_in,g)
    D_score = D_score_in; % <--- 关键：在函数内部创建 D_score 的副本
    pop=zeros(popnum,D);
    gg=cell2mat(g);
    for i=1:popnum
        mustselectnum_position=randperm(size(gg,1),1);
        mustselectnum=gg(mustselectnum_position,:);
        pop(i,mustselectnum)=1;
        D_score(mustselectnum)=0;
        canditateD=find(D_score~=0);

        min_D_score=D_score;
        gennum=round(0.9*rand*length(canditateD));
        for j=1:gennum
            canditateD=find(min_D_score~=0);
            N = numel(canditateD);
            if N < 2
                break;
            end
            variables = randperm(N,2);
            if D_score(canditateD(variables(1)))>D_score(canditateD(variables(2)))
                pop(i,canditateD(variables(1)))=1;
                min_D_score(canditateD(variables(1)))=0;
            else
                pop(i,canditateD(variables(2)))=1;
                min_D_score(canditateD(variables(2)))=0;
            end
        end
    end
end

% function [OffMask,poss_s_num] = Operator_adrpv(ParentMask,Site,rbm,allOne,other,D_score,test_adjacency,g,pv,UpdateRatio,current_subspace)
%     % 入参说明：
%     %   ParentMask - 全维度父代池（N_parent×D_full，原始父代掩码）
%     %   current_subspace - 子空间维度索引（1×D）
%     % 其他入参不变
% 
%     D = size(ParentMask, 2);
%     S = length(current_subspace);  % 子空间维度数
%     other = find(other);  % 将other逻辑数组转为索引数组（哪些列是目标列）
%     other_sub = intersect(other, current_subspace);% other_sub：current_subspace中属于other的列（索引数组，列数=目标维度数）  
%     other_nonsubspace = setdiff(other, other_sub);
%     %实际上需要用交叉变异算子Operator改动的只有other部分(other_sub部分)
% 
%     % 父代拆分（全维度）
%     Parent1Mask = ParentMask(1:floor(end/2), :);
%     Parent2Mask = ParentMask(floor(end/2)+1:end, :);
%     Site1 = Site(:,1:floor(end/2));
%     Site2 = Site(:,floor(end/2)+1:end);
% 
%     OffMaskT1 = [];
% 
%     % =======================
%     % ① RBM 分支
%     % =======================
%     if any(Site)
%         % if isempty(other_sub)
%         %     other_sub = current_subspace;
%         % end
% 
%         % ============== RBM 生成子代 ==============
%         poss = rbm.reduce(ParentMask(:,other));   % safe now
%         [OffTemp,parent1_indices] = Cross2(sum(Site),size(poss,2),poss);
%         OffTemp = rbm.recover(OffTemp);
% 
%         OffMaskT1 = false(size(OffTemp,1), D);
%         OffMaskT1(:,allOne) = true;
%         OffMaskT1(:,other) = OffTemp;
% 
%         % ParentMask_Site = ParentMask(Site,:);
%         % OffMaskT1(:,other_nonsubspace) = ParentMask_Site(parent1_indices, other_nonsubspace);
% 
%     end
% 
%     %% 非RBM父代的Cross_pv逻辑（同理改造，确保非子空间继承Parent1Mask）
%     OffMaskT2 = false(floor(sum(~Site)/2),size(ParentMask,2));
%     if ~isempty(other_sub)
% 
%         [OffMaskT2] = Cross_pv(Parent1Mask(~Site1,:), Parent2Mask(~Site2,:), other, pv, current_subspace);
%     end
%     %OffMaskT2(:, allOne) = true;
% 
%     %% 合并子代+变异（Muation_T仅子空间变异）
%     OffMask = [OffMaskT1; OffMaskT2];
%     if isempty(OffMask)
%         poss_s_num = 0;
%         return;
%     end
%     [OffMask,poss_s_num] = Muation_T(OffMask,D_score,test_adjacency,Site);
% end

