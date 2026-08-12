function [OffMask,poss_s_num] = Operator_adr(ParentMask,Site,rbm,allOne,other,D_score,test_adjacency,g,UpdateRatio,current_subspace) 
    D = size(ParentMask, 2);
    % --- 💡 关键修正：确保 current_subspace 是索引数组 ---
    if islogical(current_subspace)
        subspace_idx = find(current_subspace); 
    else
        subspace_idx = current_subspace;
    end
    
    % 确保 other 也是索引数组
    if islogical(other)
        other_idx = find(other);
    else
        other_idx = other;
    end
    
    % 1. 计算需要在子空间内交叉的列（索引）
    % intersect 会返回 other 中且属于当前子空间的索引
    other_sub = intersect(other_idx, subspace_idx); 
    
    % 2. 计算在 other 中但不在当前子空间内的列（保持父代原样）
    other_nonsubspace = setdiff(other_idx, other_sub);
    
    % 如果 other_sub 为空（即当前子空间没覆盖到任何目标维度），强制至少搜索子空间本身
    if isempty(other_sub)
        other_sub = subspace_idx;
    end

    OffMaskT1 = [];
    % =======================
    % ① RBM 分支 (Site 为 true 的个体)
    % =======================
    if any(Site)
        % 维度检查
        if isempty(rbm) || isempty(rbm.Weight)
            use_rbm = false;
        elseif size(ParentMask(Site,other_sub), 2) ~= size(rbm.Weight, 1)
            use_rbm = false;
        else
            use_rbm = true;
        end

        if use_rbm
            X = ParentMask(Site,other_sub);
            poss = rbm.reduce(X);
            [OffTemp, parent1_indices] = Cross2(sum(Site), size(poss,2), poss);
            OffTemp = rbm.recover(OffTemp);

            OffMaskT1 = false(size(OffTemp,1), D);
            OffMaskT1(:,allOne) = true;
            % 这里的 OffTemp 列数必须严格等于 length(other_sub)
            OffMaskT1(:,other_sub) = OffTemp; 
            
            ParentMask_Site = ParentMask(Site,:);
            OffMaskT1(:,other_nonsubspace) = ParentMask_Site(parent1_indices, other_nonsubspace);
        else
            % Fallback
            ParentMask_Site = ParentMask(Site,:);
            % 💡 修正点：Cross2 的第二个参数应为 length(other_sub)
            [OffTemp1, parent1_indices_ns] = Cross2(sum(Site), length(other_sub), ParentMask_Site(:,other_sub));
            
            OffMaskT1 = false(size(OffTemp1,1), D);
            OffMaskT1(:,allOne) = true;
            % 此处赋值必须对齐
            OffMaskT1(:,other_sub) = OffTemp1; 
            OffMaskT1(:,other_nonsubspace) = ParentMask_Site(parent1_indices_ns, other_nonsubspace);
        end
    end
    % =======================
    % ② 非 RBM 分支 (Site 为 false 的个体)
    % =======================
    ParentMask_nSite = ParentMask(~Site,:);
    [OffTemp2, parent1_indices_ns2] = Cross2(sum(~Site), length(other_sub), ParentMask_nSite(:,other_sub));
    
    OffMaskT2 = false(size(OffTemp2,1), D);
    OffMaskT2(:,allOne) = true;
    OffMaskT2(:,other_sub) = OffTemp2;
    OffMaskT2(:,other_nonsubspace) = ParentMask_nSite(parent1_indices_ns2, other_nonsubspace);
    
    OffMask = [OffMaskT1; OffMaskT2];
       
    [OffMask, poss_s_num] = Muation_T(OffMask, D_score, test_adjacency, Site);
end
function pop=creatpop(popnum,D,D_score,g)
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

% function [OffMask,poss_s_num] = Operator_adr(ParentMask,Site,rbm,allOne,other,D_score,test_adjacency,g,UpdateRatio,current_subspace) 
%     D = size(ParentMask, 2);
%     % --- 💡 关键修正：确保 current_subspace 是索引数组 ---
%     if islogical(current_subspace)
%         subspace_idx = find(current_subspace); 
%     else
%         subspace_idx = current_subspace;
%     end
% 
%     % 确保 other 也是索引数组
%     if islogical(other)
%         other_idx = find(other);
%     else
%         other_idx = other;
%     end
% 
%     % 1. 计算需要在子空间内交叉的列（索引）
%     % intersect 会返回 other 中且属于当前子空间的索引
%     other_sub = intersect(other_idx, subspace_idx); 
% 
%     % 2. 计算在 other 中但不在当前子空间内的列（保持父代原样）
%     other_nonsubspace = setdiff(other_idx, other_sub);
% 
%     % 如果 other_sub 为空（即当前子空间没覆盖到任何目标维度），强制至少搜索子空间本身
%     if isempty(other_sub)
%         other_sub = subspace_idx;
%     end
% 
%     OffMaskT1 = [];
%     % =======================
%     % ① RBM 分支 (Site 为 true 的个体)
%     % =======================
%     if any(Site)
%         poss = rbm.reduce(ParentMask(:,other));
%         [OffTemp, parent1_indices] = Cross2(sum(Site), size(poss,2), poss);
%         OffTemp = rbm.recover(OffTemp);
% 
%         OffMaskT1 = false(size(OffTemp,1), D);
%         OffMaskT1(:,allOne) = true;
%         % 这里的 OffTemp 列数必须严格等于 length(other_sub)
%         OffMaskT1(:,other) = OffTemp;
% 
%         % ParentMask_Site = ParentMask(Site,:);
%         % OffMaskT1(:,other_nonsubspace) = ParentMask_Site(parent1_indices, other_nonsubspace);
%     end
%     % =======================
%     % ② 非 RBM 分支 (Site 为 false 的个体)
%     % =======================
%     ParentMask_nSite = ParentMask(~Site,:);
%     [OffTemp2, parent1_indices_ns2] = Cross2(sum(~Site), length(other_sub), ParentMask_nSite(:,other_sub));
% 
%     OffMaskT2 = false(size(OffTemp2,1), D);
%     OffMaskT2(:,allOne) = true;
%     OffMaskT2(:,other_sub) = OffTemp2;
%     OffMaskT2(:,other_nonsubspace) = ParentMask_nSite(parent1_indices_ns2, other_nonsubspace);
% 
%     OffMask = [OffMaskT1; OffMaskT2];
% 
%     [OffMask, poss_s_num] = Muation_T(OffMask, D_score, test_adjacency, Site);
% end



