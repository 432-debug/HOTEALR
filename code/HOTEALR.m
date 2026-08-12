% =========================================================================
% DSSEA_TSSKEA_KEGG 独立运行版本
% 阶段 1: 用于降维的变量评分+聚类
% 阶段 2: MP-MMEA 多种群精细搜索
% =========================================================================

function [FVV, Final_Pop, POP_all, FV_all, Non_dominated_sol_MPMMEA] = Motifsubspace_Nodenew(maxFE_total, N, path_in)
    CANCER=cell(3,1);
    CANCER{1,1}='BRCA';
    CANCER{2,1}='LUSC';
    CANCER{3,1}='LUAD';

    path_in='F:\YAN1\PDENB-Copy\Motif_subspace\TSSKEA-main\Data\';
    path_out='F:\YAN1\PDENB-Copy\Motif_subspace\Motifsubspace_BRCA_Node\Motifsubspace_Node1000_';
    motif_path = 'F:\YAN1\Motif\chain_Motifcluster1node_BRCA_result';
    
    for cn=1:1
        cancer=CANCER{cn,1};
        EXP_num=30;
        Problem.maxFE_total=1000;
        Problem.N=300;
        Problem.lower=0;
        Problem.upper=1;
        for fnum = 109 :109
            if cn==3
                fnum=fnum+49;
            end
            % --- 增加总计时器 ---
            total_runtime_start = tic;
            % --- 参数配置 ---
            data=load([path_in,cancer,'_node_net\DE_BRCA_sample_result_',num2str(fnum),'.mat']);
            test_adjacency = data.subnetwork_adjacency;
            %gg=data.gg;
            kegg_file = 'KEGG_244_pathways_network.xlsx';
            [pair2path, path2pairs, pair_score] = init_pathway_mapping_node(kegg_file, data.subnetwork_genes);
    
            Problem.D = size(test_adjacency, 1); % 维度
            test_net = zeros(Problem.D, Problem.D);
            test_net(test_adjacency ~= 0) = 1;
    
            Cons = CONS(test_net);
            Cons(all(Cons == 0, 2), :) = [];
            Cnum = size(Cons, 1);
            g = cell(Cnum, 1);
            for i = 1:Cnum
                g{i} = find(Cons(i, :) == 1);
            end
            Dimension = eye(Problem.D, Problem.D);
            CV_num = Calcons(Problem.D, Cnum, g, Dimension);
            D_score = abs(CV_num - Cnum); %D_score越高节点越重要
            D_score_mean = mean(D_score);
            D_score_std  = std(D_score);

            fprintf('>>> D_score 统计 | mean = %.4f, std = %.4f, max = %.4f, min = %.4f\n', ...
                D_score_mean, D_score_std, max(D_score), min(D_score));

            motif_file = fullfile(motif_path, ...
                ['DE_BRCA_node_biomarker_', num2str(fnum), '_Motifcluster_chain.mat']);
            fprintf("== 读取 Motif 聚类文件: %s ==\n", motif_file);
            motif_data  = load(motif_file);
            motif_label = motif_data.motif_label_super;
            chain_vec   = motif_data.chain_vec;
            Non_dominated_sol=cell(30,1);
            for EXP_NUM=1:EXP_num
                Current_FE = 0;
                %% ------------------- 阶段 2: MP-MMEA 精细搜索 ---------------------------------
                fprintf('--- 阶段 2: MP-MMEA 精细搜索开始 ---\n');
                PV = tiedrank(-D_score)'; %pv越大越好，越重要
    
                [subMasks, All_Masks, initial_subspace] = motif_vssps_enhanced(motif_label, chain_vec, test_adjacency, Problem.N, D_score);
                K = length(subMasks); 
                fprintf('MP-MMEA子种群数：K=%d（与AP有效聚类数一致）\n', K);
    
                % --- 步骤3：初始化子种群参数 ---
                Masks       = cell(1,K);
                Decs        = cell(1,K);
                Populations = cell(1,K);
                GV          = cell(1,K);
                FrontNo     = cell(1,K);
                CrowdDis    = cell(1,K);
                sv          = cell(1,K);
                pv          = cell(1,K);
                Last_temp_num = 0;
                current_subspace = cell(1,K);
    
                for i = 1 : K
                    Masks{i} = subMasks{i};
                    Decs{i} = ones(size(Masks{i}));
                    N_sub_i = size(Masks{i}, 1); 
    
                    [Populations{i}, calnum_off0] = SOLUTION(Masks{i}, test_adjacency, g, D_score);
                    calnum_off0 = norm_calnum(calnum_off0, 'InitPop');
                    Current_FE = double(Current_FE) + calnum_off0;
                    
                    [Populations{i}, Decs{i}, Masks{i}, FrontNo{i}, CrowdDis{i}] = ...
                        EnvironmentalSelection(Populations{i}, ones(N_sub_i, Problem.D), Masks{i}, N_sub_i);
    
                    GV{i} = UpdateGV(zeros(1, Problem.D), Masks{i}, FrontNo{i});
                    sv{i} = zeros(1, Problem.D);
                    pv{i} = PV;
                    current_subspace{i} = zeros(1,Problem.D);
                    
                    fprintf('  子种群%d：%d个个体，掩码维度%d×%d\n', ...
                        i, N_sub_i, size(Masks{i},1), size(Masks{i},2));
                end
                % =========================================================
                % 🧹 新增：初始化后立即清理空种群 🧹
                % =========================================================
                empty_idx = [];
                for k_chk = 1:K
                    if isempty(Populations{k_chk}) || size(Masks{k_chk},1) < 2
                        empty_idx = [empty_idx, k_chk];
                    end
                end
                if ~isempty(empty_idx)
                    fprintf('初始化检测：移除 %d 个无效/空子种群。\n', length(empty_idx));
                    Populations(empty_idx) = [];
                    Masks(empty_idx) = [];
                    Decs(empty_idx) = [];
                    GV(empty_idx) = [];
                    FrontNo(empty_idx) = [];
                    CrowdDis(empty_idx) = [];
                    sv(empty_idx) = [];
                    pv(empty_idx) = [];
                    current_subspace(empty_idx) = [];
                    K = length(Populations); % 更新 K
                end
                % =========================================================
                Non_dominated_sol_MPMMEA = {};
                FE_cycle_count = 0;
                rho = 0.5;
                RHO=0.5;
                sv_global = zeros(1, Problem.D);
                NodeScore = zeros(1, Problem.D);   % 节点长期价值
                alpha_ns = 0.9;                   % 遗忘因子

                % =========================================================
                % 【关键修复】 将停滞检测变量初始化移到 while 循环之前
                % =========================================================
                Last_Min_Obj = ones(1, K) * inf; % 记录每个子种群上一代的最小目标值
                Stagnation_Count = zeros(1, K);  % 停滞计数器
                Restart_Threshold = 15;          % 设定停滞代数阈值
    
                %% ------------------- 步骤5：动态子空间下的MP-MMEA种群迭代 -------------------
                while (Current_FE < Problem.maxFE_total)
                    FE_cycle_count = FE_cycle_count + 1;
                    delta= Current_FE/Problem.maxFE_total;
                        % 1. 收集所有子种群的 F1 个体 Mask
                        F1_Parents_Masks = [];
                        K = length(Populations); 
                        for j = 1:K
                            if ~isempty(FrontNo{j}) && ~isempty(Masks{j})
                                F1_idx = (FrontNo{j} == 1);
                                F1_Parents_Masks = [F1_Parents_Masks; Masks{j}(F1_idx, :)];
                            end
                        end
                        if isempty(F1_Parents_Masks)
                            Parent_Masks_for_KSS = OffMask; 
                        else
                            Parent_Masks_for_KSS = F1_Parents_Masks;
                        end                
    
                    % 子种群排序 
                    [~,rank] = sort(SubPopRank(Populations)); 
                    for i = 1 : K
                        sub_idx = rank(i);
                        N_sub_i = length(Populations{sub_idx});

                        % 【修复代码】：如果子种群个体数少于2，无法进行锦标赛选择，直接跳过
                        % 这种种群会在稍后的“种群崩溃合并逻辑”中被清理掉
                        if N_sub_i < 2
                            % fprintf('  DEBUG: 跳过子种群 %d (N=%d)\n', sub_idx, N_sub_i);
                            continue;
                        end

                        GV{sub_idx} = UpdateGV(GV{sub_idx},Masks{sub_idx},FrontNo{sub_idx});
                        First_Masks{sub_idx}=Masks{sub_idx}(FrontNo{sub_idx}==1,:);
                        [temp_num,~]=size(First_Masks{sub_idx});
                        temp_vote=sum(First_Masks{sub_idx},1);
                        sv{sub_idx}(1,:)=(Last_temp_num/(Last_temp_num+temp_num))*sv{sub_idx}(1,:)+(temp_num/(Last_temp_num+temp_num))*(temp_vote/temp_num);
                        Last_temp_num=temp_num;
                        eta_i = temp_num / N_sub_i;
                        pv{sub_idx} = pv{sub_idx} + 0.5 * (eta_i + delta) .* (1 - sv{sub_idx}) .* pv{sub_idx};
                        % 归一化pv{i}（保持与文献一致的更新后处理）
                        pv{sub_idx} = pv{sub_idx} / max(pv{sub_idx} + eps);  % 避免除以零
                        % ===== 新增：全局 sv_global（只累积 F1）=====
                        % 收集所有子种群 F1（一次 per FE_cycle）
                        F1_all = [];
                        for kk = 1:K
                            if isempty(Masks{kk}) || isempty(FrontNo{kk}), continue; end
                            idx_f1 = (FrontNo{kk} == 1);
                            if any(idx_f1)
                                F1_all = [F1_all; Masks{kk}(idx_f1,:)];
                            end
                        end

                        if ~isempty(F1_all)
                            vote_global = mean(F1_all,1);
                            alpha = 0.2;   % 推荐 0.1~0.3，重启用不宜太大
                            sv_global = (1-alpha)*sv_global + alpha*vote_global;
                        end
                        
                        current_subspace{sub_idx} = dynamic_subspace_pareto_neighbors(Masks{sub_idx}, FrontNo{sub_idx}, test_adjacency);
                        
                        Mating = TournamentSelection(2, N_sub_i, FrontNo{sub_idx}, -CrowdDis{sub_idx});
                        MatingMasks = Masks{sub_idx}(Mating,:);
    
                        % Operator
                        if (delta/0.618) < 0.618 %% 第一阶段pv
                            Site=false(1,N_sub_i);
                            [rbm,allOne,other] = deal([]);
                            [OffMask,poss_s_num] = Operator_adrpv(MatingMasks,Site,rbm,allOne,other,...
                                D_score, test_adjacency, g, pv{sub_idx}, rho, current_subspace{sub_idx});                        
                        else %% 第二阶段pv+sv+cv
                            if rho<0.5
                                index=FrontNo{sub_idx}<ceil(max(FrontNo{sub_idx})/2);
                            else
                                index=FrontNo{sub_idx}==1;
                            end
                            if size(Masks{sub_idx},1) < 10
                                rbm = []; allOne = []; other = [];
                                Site = false(size(Masks{sub_idx},1),1); 
                            else
                                [rbm,allOne,other,Site] = trainRBM(Masks{sub_idx}(index,:),rho,size(Masks{sub_idx},1),current_subspace{sub_idx});
                            end
    
                            if rand < 1/3
                                [OffMask,poss_s_num] = Operator_adrpv(MatingMasks,Site,rbm,allOne,other,...
                                    D_score,test_adjacency,g,pv{sub_idx},rho,current_subspace{sub_idx});
                            elseif rand > 2/3
                                [OffMask,poss_s_num] = Operator_adrsv(MatingMasks,Site,rbm,allOne,other, ...
                                    D_score,test_adjacency,g,sv{sub_idx},rho,current_subspace{sub_idx});
                            else
                                [OffMask,poss_s_num] = Operator_adr(MatingMasks,Site,rbm,allOne,other, ...
                                    D_score,test_adjacency,g,rho,current_subspace{sub_idx});
                            end
                        end
    
                        KSS = Extract_KSS_Cooccur2(Parent_Masks_for_KSS, test_adjacency, D_score);
                        if ~exist('KSS', 'var') || isempty(KSS), KSS = {}; end
                        New_OffMask = Operator_KSS_Guided2(Parent_Masks_for_KSS,KSS,D_score,test_adjacency);                    
                        OffMask = [New_OffMask; OffMask];
    
                        OffDec = ones(size(OffMask));
                        [Offspring, calnum_off1] = SOLUTION(OffMask, test_adjacency, g, D_score);
                        calnum_off1 = norm_calnum(calnum_off1, 'Offspring');
                        Current_FE = double(Current_FE) + calnum_off1;
    
                        % Populations{sub_idx} = [Populations{sub_idx}, Offspring];
                        % Decs{sub_idx} = [Decs{sub_idx}; OffDec];
                        % Masks{sub_idx} = [Masks{sub_idx}; OffMask];
    
                        if i > 1 
                            R = zeros(1,Problem.D);
                            for j = 1 : i-1
                                best_sub_idx = rank(j);
                                R = R + GV{best_sub_idx}; 
                            end
                            R(R>0) = 1;
                            dis = sum(repmat(R,length(Populations{rank(i)}),1)&Masks{rank(i)},2);
                            [Populations{sub_idx},Decs{sub_idx},Masks{sub_idx},FrontNo{sub_idx},CrowdDis{sub_idx},sRatio] = EnvironmentalSelection_adr( ...
                                [Populations{sub_idx}, Offspring],[Decs{sub_idx}; OffDec],[Masks{sub_idx}; OffMask],N_sub_i,length(Populations{sub_idx}),poss_s_num); 
                        else
                             [Populations{sub_idx},Decs{sub_idx},Masks{sub_idx},FrontNo{sub_idx},CrowdDis{sub_idx},sRatio] = EnvironmentalSelection_adr( ...
                                 [Populations{sub_idx}, Offspring],[Decs{sub_idx}; OffDec],[Masks{sub_idx}; OffMask],N_sub_i,length(Populations{sub_idx}),poss_s_num);
                        end
                        rho = (rho+sRatio)/2;
                        RHO=[RHO;rho];
                        if rand < 0.5 
                            CV_i = extract(Populations{i});
                            [FrontNo_i,~] = NDSort(CV_i, size(Populations{i}, 2));
                            first_front_idx = find(FrontNo_i == 1);
            
                            candidate_masks = Masks{i}(first_front_idx, :);
                            LS_masks = local_search_pathway_node(candidate_masks, data, pair2path, path2pairs, pair_score);
                            [LS_pop, calnum_LS] = SOLUTION(LS_masks, test_adjacency, g, D_score);
                            calnum_LS = norm_calnum(calnum_LS, 'LocalSearch');
                            Current_FE = double(Current_FE) + calnum_LS;
                            
                            combinedPop = [Populations{i}, LS_pop];
                            combinedMasks = [Masks{i}; logical(LS_masks)];
                            combinedDecs = [Decs{i}; ones(size(LS_masks))];
                            [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection...
                                (combinedPop, combinedDecs, combinedMasks, N_sub_i);                     
                        end
                    end % for i = 1 : K
    
                    % =========================================================
                    % 📢 种群个体数崩溃合并逻辑 (N_sub < 20)
                    % =========================================================
                    N_min_threshold = 20; 
                    if K > 1 
                        current_pop_sizes = cellfun(@length, Populations);
                        target_indices = find(current_pop_sizes < N_min_threshold);
                        if ~isempty(target_indices)
                            fprintf('⚠️ 检测到 %d 个子种群个体数低于 %d，执行合并。\n', length(target_indices), N_min_threshold);
                            for target_idx_k = 1:length(target_indices)
                                i_target = target_indices(target_idx_k); 
                                if i_target > K || isempty(Populations{i_target}), continue; end
                                valid_indices = 1:K;
                                valid_indices(i_target) = [];
                                if isempty(valid_indices), continue; end
                                [~, min_size_rel_idx] = min(current_pop_sizes(valid_indices));
                                j_partner = valid_indices(min_size_rel_idx); 
                                
                                fprintf(' 合并子种群 %d (%d 个个体) 到子种群 %d (%d 个个体)。\n', ...
                                    i_target, length(Populations{i_target}), ...
                                    j_partner, length(Populations{j_partner}));
                    
                                Populations{j_partner} = [Populations{j_partner}, Populations{i_target}];
                                Decs{j_partner}        = [Decs{j_partner}; Decs{i_target}];
                                Masks{j_partner}       = [Masks{j_partner}; Masks{i_target}];
                                if i_target < j_partner, j_partner = j_partner - 1; end
                                
                                % 删除被合并种群数据
                                Populations(i_target) = []; Decs(i_target) = []; Masks(i_target) = []; 
                                GV(i_target) = []; sv(i_target) = []; pv(i_target) = []; FrontNo(i_target) = []; CrowdDis(i_target) = [];
                                current_subspace(i_target) = [];
                                
                                % 【同步删除停滞检测变量】
                                Last_Min_Obj(i_target) = [];
                                Stagnation_Count(i_target) = [];
                                
                                K = K - 1;
                                GV{j_partner} = UpdateGV(zeros(1, Problem.D), Masks{j_partner}, FrontNo{j_partner});
                                current_pop_sizes = cellfun(@length, Populations);
                                [~, rank] = sort(SubPopRank(Populations));
                                break; 
                            end 
                        end 
                    end 
                    
                    % =========================================================
                    % 6. 合并-划分操作 (简化: 每 50 次循环)
                    % =========================================================
                    if mod(FE_cycle_count, 30)==0
                        [~,best_rank] = SubPopRank(Populations);
                        divisionFlag = all(best_rank==1);
                        [ss,index]   = SubPopSimility(Populations,Masks);
                    
                        if ss > 0.5 && K > 1 % 合并操作
                            i = index(1); j = index(2); 
                            [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection( ...
                                [Populations{i},Populations{j}],[Decs{i};Decs{j}],[Masks{i};Masks{j}],floor(Problem.N/K));
                
                            Populations(j) = []; Decs(j) = []; Masks(j) = []; GV(j) = []; sv(j) = []; pv(j) = [];
                            FrontNo(j) = []; CrowdDis(j) = []; current_subspace(j) = []; 
                            
                            % 【同步删除停滞检测变量】
                            Last_Min_Obj(j) = [];
                            Stagnation_Count(j) = [];
                            
                            K = K-1;
                            N_sub_i = floor(Problem.N/K);
                            fprintf('--- 合并操作: K=%d ---\n', K);
                        
                        elseif divisionFlag == 1 % 划分操作
                            % 新增：如果当前种群太小，不要强行划分，否则下一轮又被合并了
                            if N_sub_i < 40
                                continue;
                            end
                            N_sub_new = floor(Problem.N/(K+1));
                            for i_div = 1 : K
                                [Populations{i_div},Decs{i_div},Masks{i_div},FrontNo{i_div},CrowdDis{i_div}] = EnvironmentalSelection( ...
                                    Populations{i_div},Decs{i_div},Masks{i_div},N_sub_new);
                            end
                            K = K + 1;
                            Dec_new  = ones(N_sub_new, Problem.D);            
                            F = zeros(1,Problem.D);
                            for i_div = 1: K-1, F = F + GV{i_div}; end            
                            Mask_new = zeros(N_sub_new, Problem.D);
                            for i_new = 1 : N_sub_new
                                candidate_vars = find(F == min(F)); 
                                if isempty(candidate_vars), candidate_vars = 1:Problem.D; end 
                                num_select = randi([1, min(10, length(candidate_vars))]); 
                                idx_select = randperm(length(candidate_vars), num_select);
                                Mask_new(i_new, candidate_vars(idx_select)) = 1;
                                gg=cell2mat(g);
                                mustselectnum=gg(randperm(size(gg,1),1),:);
                                Mask_new(i_new,mustselectnum(mustselectnum~=0))=1;
                                Mask_new(i_new,:) = logical(Mask_new(i_new,:));
                            end
                    
                            [Populations{K}, calnum_new] = SOLUTION(Mask_new, test_adjacency, g, D_score);
                            calnum_new = norm_calnum(calnum_new, 'Division');
                            Current_FE = double(Current_FE) + calnum_new;
    
                            Masks{K} = Mask_new; Decs{K} = Dec_new; GV{K} = zeros(1,Problem.D); sv{K} = zeros(1,Problem.D); pv{K} = zeros(1,Problem.D);
                            current_subspace{K} = zeros(1, Problem.D); 
                
                            [Populations{K},Decs{K},Masks{K},FrontNo{K},CrowdDis{K}] = EnvironmentalSelection(Populations{K},Decs{K},Masks{K},length(Populations{K}));
                            GV{K} = UpdateGV(zeros(1,Problem.D),Masks{K},FrontNo{K});
                            
                            % 【同步增加停滞检测变量】
                            Last_Min_Obj(K) = inf;
                            Stagnation_Count(K) = 0;
                            
                            N_sub = floor(Problem.N/K);
                            fprintf('--- 划分操作: K=%d ---\n', K);
                        end
                    end
                    % 7. 汇总当前所有非支配解
                    Population_all = [Populations{:}];
                    FV_all = extract(Population_all);
                    [FrontNo_all, ~] = NDSort(FV_all, size(FV_all, 1));
                    Pop_final = extract_d(Population_all(FrontNo_all == 1));
                    FV_all = FV_all(FrontNo_all == 1, :);
    
                    if isempty(Current_FE) || ~isnumeric(Current_FE)
                        cf = -1;
                    else
                        cf = round(Current_FE);
                    end
                    % ======= 计算NodeScore=======
                    for k = 1:length(Populations)
                        if isempty(Populations{k}), continue; end
                        F1_idx = (FrontNo{k} == 1);
                        if any(F1_idx)
                            F1_Masks = Masks{k}(F1_idx,:);
                            NodeScore = alpha_ns * NodeScore + (1-alpha_ns) * mean(F1_Masks,1);
                        end
                    end
                    % =====================================================================
                    % 🛑 防早熟机制：停滞重启 + 多样性注入 (Stagnation Restart) 🛑
                    % =====================================================================
                    for k_chk = 1 : length(Populations) 
                        if isempty(Populations{k_chk}), continue; end
    
                        % 【修复点】: 跳过过小的种群，避免索引错误
                        N_curr = length(Populations{k_chk});
                        if N_curr < 10, continue; end 
    
                        Current_Objs = extract(Populations{k_chk});
                        Current_Min_Obj = min(sum(Current_Objs, 2));
                        
                        % 防止 Last_Min_Obj 维度不匹配（防御性编程）
                        if k_chk > length(Last_Min_Obj)
                            Last_Min_Obj(k_chk) = inf;
                            Stagnation_Count(k_chk) = 0;
                        end
    
                        if abs(Current_Min_Obj - Last_Min_Obj(k_chk)) < 1e-4
                            Stagnation_Count(k_chk) = Stagnation_Count(k_chk) + 1;
                        else
                            Stagnation_Count(k_chk) = 0; 
                            Last_Min_Obj(k_chk) = Current_Min_Obj;
                        end
    
                        if Stagnation_Count(k_chk) >= Restart_Threshold
                            fprintf('  >> [重启] 子种群 %d 已停滞 %d 代，执行多样性注入...\n', k_chk, Stagnation_Count(k_chk));
                            
                            Keep_Num = floor(N_curr * 0.5);
                            
                            % 【修复点】: 使用 FrontNo 和 CrowdDis 排序，替代 SubPopRank
                            % 确保为列向量
                            F_col = reshape(FrontNo{k_chk}, [], 1);
                            C_col = reshape(CrowdDis{k_chk}, [], 1);
                            % 优先 FrontNo (升序), 其次 CrowdDis (降序)
                            [~, Rank_Local] = sortrows([F_col, -C_col]);
                            
                            % 截取精英
                            Keep_Ind = Rank_Local(1:Keep_Num);
                            Pop_Keep = Populations{k_chk}(Keep_Ind);
    
                            External_Elites = [];
                            for other_k = 1 : length(Populations)
                                if other_k ~= k_chk && ~isempty(Populations{other_k})
                                    % 确保维度安全
                                    if length(FrontNo{other_k}) == length(Populations{other_k})
                                        other_F1_Idx = (FrontNo{other_k} == 1);
                                        if any(other_F1_Idx)
                                            External_Elites = [External_Elites, Populations{other_k}(other_F1_Idx)];
                                        end
                                    end
                                end
                            end
                            Inject_Num = N_curr - Keep_Num;
                            Pop_Inject = {};
    
                            if ~isempty(External_Elites)
                                Num_From_External = floor(Inject_Num * 0.5);
                                Rand_Idx = randperm(length(External_Elites), min(Num_From_External, length(External_Elites)));
                                Pop_Inject = [Pop_Inject, External_Elites(Rand_Idx)];
                            end
                            Num_Random = Inject_Num - length(Pop_Inject);
                            if Num_Random > 0
                                Rand_Mask_All = [];

                                Num_Walk  = floor(Num_Random * 0);
                                Num_VSSPS = Num_Random - Num_Walk;
                                if Num_Walk > 0
                                    % Walk_Mask = rand(Num_Walk, Problem.D) < 0.05;
                                    p_start = NodeScore + 1e-6;
                                    p_start = p_start / sum(p_start);
                                    Walk_Mask = false(Num_Walk, Problem.D);
                                    for r = 1:Num_Walk
                                        curr = randsample(Problem.D, 1, true, p_start);
                                        Walk_Mask(r, curr) = true;
                                        % 简单扩散1-2步
                                        steps = randi([2,5]);
                                        for s = 1:steps
                                            act = find(Walk_Mask(r,:));
                                            [~, nbs] = find(test_adjacency(act, :));
                                            nbs = unique(nbs);
                                            if ~isempty(nbs)
                                                Walk_Mask(r, nbs(randi(length(nbs)))) = true;
                                            end
                                        end
                                    end
                                    [Pop_Random, calnum_rnd] = SOLUTION(Walk_Mask, test_adjacency, g, D_score);
                                    calnum_rnd = norm_calnum(calnum_rnd, 'Restart');
                                    Current_FE = Current_FE + calnum_rnd;
                                    Pop_Inject = [Pop_Inject, Pop_Random];
                                end

                                % ========== 策略 B（替换）：DSMLBEM 层级进化多样性搜索 ==========
                                if Num_VSSPS > 0
                                    % -------- 全局个体池 --------
                                    All_Pop   = [];
                                    All_Masks = [];
                                    All_Obj   = [];
                                    for kk = 1:length(Populations)
                                        if isempty(Populations{kk}), continue; end
                                        All_Pop = [All_Pop, Populations{kk}(:)'];
                                        All_Masks = [All_Masks; Masks{kk}];
                                        All_Obj   = [All_Obj; extract(Populations{kk})];
                                    end

                                    WarningScore = All_Obj(:,2);   % 你已有的假设
                                    [~, idx_sort] = sort(WarningScore, 'ascend');

                                    Np_all = size(All_Masks,1);
                                    layer_size = floor(Np_all / 8);

                                    Layers = cell(8,1);
                                    for l = 1:8
                                        if l < 8
                                            Layers{l} = idx_sort((l-1)*layer_size+1 : l*layer_size);
                                        else
                                            Layers{l} = idx_sort((l-1)*layer_size+1 : end);
                                        end
                                    end

                                    % --- 3. 定义协同进化层级对 ---
                                    CoopPairs = {
                                        1, 2;   % L1 × L1
                                        3, 4;   % L1 × L2
                                        %3, 4;   % L1 × L6
                                        };
                                    pid = randi(size(CoopPairs,1));
                                    La = CoopPairs{pid,1};
                                    Lb = CoopPairs{pid,2};

                                    New_Masks = zeros(Num_VSSPS, Problem.D);

                                    % --- 4. 生成 DSMLBEM 子代 ---
                                    CandIdx = unique([Layers{La}; Layers{Lb}]);
                                    
                                    CandPop   = All_Pop(CandIdx);
                                    CandMasks = All_Masks(CandIdx,:);
                                    CandObj   = All_Obj(CandIdx,:);
                                    CandDecs = ones(size(CandMasks));
                                    PopObjs_cand = CandObj;   
                                    N_cand = size(PopObjs_cand,1);

                                    % 非支配排序
                                    [FrontNo_cand, ~] = NDSort(PopObjs_cand, N_cand);

                                    % 拥挤距离（注意：用 Mask 或 PopObjs 都可以）
                                    CrowdDis_cand = CrowdingDistance(CandMasks, FrontNo_cand);

                                    % Mating = TournamentSelection(2, size(CandMasks,1), FrontNo_cand, -CrowdDis_cand);
                                    % MatingMasks_restart = CandMasks(Mating,:);
                                    Mating = TournamentSelection_hamming(size(CandMasks,1),size(CandMasks,1), CandMasks, FrontNo_cand, -CrowdDis_cand);
                                    MatingMasks_restart = CandMasks(Mating,:);
                                    % -------- 4. 调用“主循环算子” --------
                                    % 仍然使用你原来的阶段控制逻辑
                                    Site = false(1, size(CandMasks,1));
                                    rbm  = []; allOne = []; other = [];
                                    if (delta/0.618) < 0
                                        % 第一阶段：pv 主导（稳）
                                        [OffMask, ~] = Operator_adrpv( ...
                                            MatingMasks_restart, Site, rbm, allOne, other, ...
                                            D_score, test_adjacency, g, PV, rho, ones(1, Problem.D));
                                    else
                                        % 第二阶段：pv / sv / adr 混合
                                        if rand < 1/3
                                            [OffMask, ~] = Operator_adrpv( ...
                                                MatingMasks_restart, Site, rbm, allOne, other, ...
                                                D_score, test_adjacency, g, PV, rho, ones(1, Problem.D));
                                        elseif rand > 2/3
                                            [OffMask, ~] = Operator_adrsv( ...
                                                MatingMasks_restart, Site, rbm, allOne, other, ...
                                                D_score, test_adjacency, g, sv_global, rho, ones(1, Problem.D));
                                        else
                                            [OffMask, ~] = Operator_adr( ...
                                                MatingMasks_restart, Site, rbm, allOne, other, ...
                                                D_score, test_adjacency, g, rho, ones(1, Problem.D));
                                        end
                                    end

                                    New_Masks = OffMask;
                                end

                                Rand_Mask_All = [Rand_Mask_All; New_Masks];
                                
                                % ========== 注入 ==========
                                if ~isempty(Rand_Mask_All)
                                    [Pop_Random, calnum_rnd] = SOLUTION(Rand_Mask_All, test_adjacency, g, D_score);
                                    calnum_rnd = norm_calnum(calnum_rnd, 'Restart');
                                    Current_FE = Current_FE + calnum_rnd;
                                    Pop_Inject = [Pop_Inject, Pop_Random];
                                end
                            end
                            % =========================================================
                            % 策略 C：非支配解低相似（Jaccard）并集生成
                            % 不设数量上限，交由环境选择控制
                            % =========================================================
                            % F1_idx = (FrontNo{k_chk} == 1);
                            F1_Masks = F1_Parents_Masks;

                            tau_jaccard = 0.4;   % 建议 0.3~0.4

                            if size(F1_Masks,1) >= 2
                                New_Masks_C = [];

                                % 两两组合
                                for i = 1:size(F1_Masks,1)-1
                                    A = F1_Masks(i,:);
                                    for j = i+1:size(F1_Masks,1)
                                        B = F1_Masks(j,:);

                                        inter = sum(A & B);
                                        uni   = sum(A | B);
                                        if uni == 0
                                            continue;
                                        end

                                        jac = inter / uni;

                                        if jac < tau_jaccard
                                            Union = A | B;
                                            keep_prob = NodeScore ./ (max(NodeScore)+eps);
                                            Mask_C = Union & (rand(1,Problem.D) < keep_prob);
                                            New_Masks_C = [New_Masks_C; Mask_C];
                                        end
                                    end
                                end

                                % --- 去重，防止完全相同的并集 ---
                                if ~isempty(New_Masks_C)
                                    New_Masks_C = unique(New_Masks_C, 'rows', 'stable');

                                    % --- 生成解并注入 ---
                                    [Pop_C, calnum_c] = SOLUTION(New_Masks_C, test_adjacency, g, D_score);
                                    calnum_c = norm_calnum(calnum_c, 'Restart');
                                    Current_FE = Current_FE + calnum_c;

                                    Pop_Inject = [Pop_Inject, Pop_C];
                                end
                            end

                            % 重组
                            % =================================================
                            % 🛡️ 稳健重组策略：精英保留 + 末位淘汰 PK 🛡️
                            % =================================================

                            % 1. 获取原种群的后半部分 (Bottom 50%)
                            % 注意：Rank_Local 是之前计算好的排序索引 [Top_50, Bottom_50]
                            Bottom_Ind = Rank_Local(Keep_Num+1 : end);
                            Pop_Bottom_Old = Populations{k_chk}(Bottom_Ind);

                            % 2. 让 [注入个体] 与 [原 Bottom 个体] 进行 PK
                            %    目标是从中选出 Inject_Num 个最好的
                            Competition_Pool = [Pop_Bottom_Old, Pop_Inject];
                            Masks_Comp = logical(extract_d(Competition_Pool));
                            Decs_Comp  = ones(size(Masks_Comp));

                            % 需要选出的数量 = 原种群被移除的数量
                            Num_To_Select = length(Bottom_Ind);

                            [Pop_Bottom_New, ~, ~, ~, ~] = ...
                                EnvironmentalSelection(Competition_Pool, Decs_Comp, Masks_Comp, Num_To_Select);

                            % 3. 最终合并：[原精英] + [PK胜出者]
                            Populations{k_chk} = [Pop_Keep, Pop_Bottom_New];

                            % 4. 更新辅助变量 (Mask, Dec, FrontNo, CrowdDis)
                            Masks{k_chk} = logical(extract_d(Populations{k_chk}));
                            Decs{k_chk}  = ones(size(Masks{k_chk}));
                            [Populations{k_chk}, Decs{k_chk}, Masks{k_chk}, FrontNo{k_chk}, CrowdDis{k_chk}] = ...
                                EnvironmentalSelection(Populations{k_chk}, Decs{k_chk}, Masks{k_chk}, length(Populations{k_chk}));

                            % 5. 更新全局变量和计数器
                            GV{k_chk} = UpdateGV(zeros(1, Problem.D), Masks{k_chk}, FrontNo{k_chk});
                            Stagnation_Count(k_chk) = 0;
                            Last_Min_Obj(k_chk) = min(sum(extract(Populations{k_chk}), 2));
                        end
                    end
                    % =====================================================================
    
                    % fprintf('FE: %d/%d, K: %d, NDS: %d\n', cf, Problem.maxFE_total, K, size(Pop_final, 1));
                    
                    if Current_FE >= Problem.maxFE_total
                        fprintf('--- 评估次数达到上限 ---\n');
                        break;
                    end
                end
                Non_dominated_sol{EXP_NUM}=Pop_final;
            end
            Non_dominated_sol1=cell2mat(Non_dominated_sol);
            [pop,~,~]=unique(Non_dominated_sol1,'rows');
            [functionvalue,~] = Calfunctionvalue(pop,test_adjacency);
            [FrontNo,~] = NDSort(functionvalue,size(functionvalue,1));
            POP=pop(FrontNo==1,:);
            FV=functionvalue(FrontNo==1,:);
            [FVV,mod_position]=sortrows(FV);
            Final_Pop=POP(mod_position,:);
            
            %% 统计多模态出现的次数c
            FVV(:,2)=-FVV(:,2);
            A=tabulate(FVV(:,1));
            B= A(:,2)~=0;
            c=A(B,:);
            %% 储存数据 
            total_runtime = toc(total_runtime_start);
            fprintf('\n=== 性能统计 ===\n');
            fprintf('总运行时间: %.4f 秒\n', total_runtime);
            disp('非支配解目标值 (FVV):');
            disp(FVV);
            filename = strcat(path_out, cancer,'_result\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_AP-DSSEA_NDSol.mat');
            save(filename, 'FVV');
            filename2 = strcat(path_out, cancer,'_result\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_AP-DSSEA_POP.mat');
            save(filename2, 'Final_Pop');
            filename3=strcat(path_out,cancer,'_result\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_AP-DSSEA_boxchart.mat');
            save(filename3,'Non_dominated_sol');
            filename4=strcat(path_out,cancer,'_result\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_AP-DSSEA_Num.mat');
            save(filename4,'c');
            Runtime_sec = total_runtime;
            filename_time = strcat(path_out, cancer,'_result\DE_',cancer,'_edge_biomarker_',num2str(fnum),'_AP-DSSEA_Time.mat');
            save(filename_time, 'Runtime_sec');
            fprintf('运行时间已保存至: %s\n', filename_time);
        end
    end
end
% --- 辅助函数保持不变 ---
function [pair2path, path2pairs, pair_score] = init_pathway_mapping_node(kegg_file, subnetwork_genes)
    % 输入：
    %   kegg_file：KEGG通路Excel文件路径
    %   subnetwork_genes：基因序号→基因名称的映射（cell数组）
    % 输出：
    %   pair2path：N×1 cell，每个节点（基因）对应的通路列表
    %   path2pairs：结构体，通路→节点索引映射
    %   pair_score：N×1向量，节点的通路共现得分

    % 1. 读取KEGG数据
    tbl = readtable(kegg_file, 'FileType', 'spreadsheet');
    if size(tbl,1) < 2
        error('KEGG数据为空或格式错误');
    end
    kegg_data = tbl{:,:};

    % 2. 确定节点数量
    n_nodes = size(subnetwork_genes, 1); 
    if n_nodes == 0
        error('subnetwork_genes为空');
    end

    % 3. 初始化 pair2path (即 node2path)
    pair2path = cell(n_nodes, 1);  
    
    % 4. 构建基因名称→序号的反向映射
    gene2idx = containers.Map('KeyType','char', 'ValueType','double');
    for i = 1:n_nodes
        gene_name = subnetwork_genes{i};
        if ~isempty(gene_name) && ~isKey(gene2idx, gene_name)
            gene2idx(gene_name) = i; 
        end
    end

    % 5. 初始化通路→节点映射结构体
    path2pairs = struct();

    % 6. 遍历KEGG数据
    for i = 1:size(kegg_data,1)
        % 提取当前行的基因A、基因B、通路名称
        g1_name = kegg_data{i,1};
        g2_name = kegg_data{i,2};
        raw_pathway = kegg_data{i,3};

        % === [关键修复] 增强的名称清洗逻辑 ===
        % 使用正则表达式将所有"非字母、非数字"的字符替换为下划线
        % 这能处理逗号(,)、冒号(:)、撇号(')、空格等所有非法字符
        pathway = regexprep(raw_pathway, '[^a-zA-Z0-9]', '_');
        
        % 将连续的下划线合并为一个 (例如 "A, B" 会变成 "A__B"，需修正为 "A_B")
        pathway = regexprep(pathway, '_+', '_');
        
        % 去除首尾可能产生的下划线
        pathway = regexprep(pathway, '^_+|_+$', '');
        
        % 检查首字符是否为数字（字段名不能以数字开头）
        if ~isempty(pathway) && isstrprop(pathway(1), 'digit')
            pathway = ['P_', pathway];
        end
        % ===================================

        % 查找基因是否存在
        indices_to_update = [];
        if isKey(gene2idx, g1_name)
            indices_to_update = [indices_to_update; gene2idx(g1_name)];
        end
        if isKey(gene2idx, g2_name)
            indices_to_update = [indices_to_update; gene2idx(g2_name)];
        end
        indices_to_update = unique(indices_to_update);

        if isempty(indices_to_update)
            continue;
        end

        % 更新映射
        for k = 1:length(indices_to_update)
            idx = indices_to_update(k);
            
            % 更新 pair2path
            if isempty(pair2path{idx})
                pair2path{idx} = {pathway}; 
            else
                if ~ismember(pathway, pair2path{idx})
                    pair2path{idx} = [pair2path{idx}; {pathway}]; 
                end
            end
            
            % 更新 path2pairs
            if isfield(path2pairs, pathway)
                path2pairs.(pathway) = unique([path2pairs.(pathway); idx]);
            else
                path2pairs.(pathway) = idx;
            end
        end
    end

    % 7. 计算得分
    pair_score = zeros(n_nodes, 1);
    for p = 1:n_nodes
        unique_paths = unique(pair2path{p});
        pair_score(p) = length(unique_paths);
    end

    % 归一化
    max_score = max(pair_score);
    if max_score > 0
        pair_score = pair_score / max_score; 
    end
end
%% ===============================================================
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

function A_adjacent=CONS(test_Net)
    [z1,z2]=find(triu(test_Net)~=0);
    z=[z1,z2];
    NNN=length(test_Net);
    N2=size(z,1);
    A_adjacent=zeros(N2,NNN);
    for i=1:N2
        A_adjacent(i,z(i,1))=1;
        A_adjacent(i,z(i,2))=1;
    end
end

function CVvalue = Calcons(popnum,Cnum,g,pop)
    cv=zeros(Cnum,1);
    CVvalue=zeros(popnum,1);
    for i=1:popnum
        ind=pop(i,:);
        for j=1:Cnum
            cv(j)=max(1-sum(ind(g{j})),0);
        end
        CVvalue(i)=sum(cv);
    end
end

function [Population,calnum]  = SOLUTION(Mask,test_adjacency,g,D_score)
% 如果 Mask 为空，返回空 population 且 calnum=0（防止上游累加变成 []）
if isempty(Mask) || size(Mask,1) == 0
    Population = {}; 
    calnum = 0;
    return;
end
[functionvalue,calnum1] = Calfunctionvalue(Mask,test_adjacency);
[Mask,tt]=FIX2(Mask,functionvalue,g,size(test_adjacency,2),D_score,test_adjacency);
[functionvalue,calnum2]=Calfunctionvalue_afterfix(Mask,test_adjacency,tt,functionvalue);
calnum=calnum1+calnum2;
Population=cell(1,size(Mask,1));%创建一个元胞数组 Population，其长度等于 Mask 的行数
for i=1:size(Population,2)
    Population{1,i}.dec=Mask(i,:);
    Population{1,i}.obj=functionvalue(i,:);
    Population{1,i}.con=0;
    Population{1,i}.add=[];
end
end
function [Functionvalue, calnum]= Calfunctionvalue(pop,test_adjacency)
Functionvalue=zeros(size(pop,1),2);
functionvalue=zeros(size(pop,1),4);
for i=1:size(pop,1)
    Functionvalue(i,1)=sum(pop(i,:));

    a=find(pop(i,:)==1);

    matrix=test_adjacency(a,:);

    inmatrix=matrix(:,a);

    genin=nonzeros(inmatrix);

    functionvalue(i,1)=abs(mean(genin));
    functionvalue(i,2)=std(genin);
    
    gen=nonzeros(matrix);
    genout=setdiff(gen,genin);
%    if isempty(genout)
%        functionvalue(i,3) = 1e-3; % 没有外部连接时，mean(genout) 设为 1e-3
%    else
%        functionvalue(i,3) = abs(mean(genout));
%    end
    functionvalue(i,3)=abs(mean(genout));
    functionvalue(i,4)=functionvalue(i,1)*functionvalue(i,2)/functionvalue(i,3);
    
end
calnum = size(pop,1);
Functionvalue(:,2)= -round(functionvalue(:,4)*100)/100;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,1)=size(test_adjacency,2);
Functionvalue(tt,2)=0;
end
function [Functionvalue,calnum] = Calfunctionvalue_afterfix(pop,test_adjacency,tt,functionvalue)

for i=1:size(tt,1)
    functionvalue(tt(i),:)=Calfunctionvalue(pop(tt(i),:),test_adjacency);
end
calnum=i;
Functionvalue=functionvalue;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,2)=0;
ttt= Functionvalue(:,2)==0;
Functionvalue(ttt,1)=size(test_adjacency,2);
end
function index = TournamentSelection_hamming(K, N, Population, FrontNo, SpCrowdDis)
% 带汉明距离辅助的锦标赛选择
% 输入：
%   K, N          - 锦标赛参数 (通常 K=2)
%   Population    - 二进制矩阵 (PopSize x D)
%   FrontNo       - 非支配排序等级 (PopSize x 1)
%   SpCrowdDis    - 拥挤距离 (PopSize x 1)
% 输出：
%   index         - 被选中的个体索引 (N x 1)

    PopSize = size(Population,1);
    if PopSize < 2
        index = ones(1, N);
        return;
    end
    % 计算汉明距离矩阵
    hamming_dist = pdist2(Population, Population, 'hamming');
    [~, site2] = sort(hamming_dist, 2);
    % 为每个选择事件生成两名候选
    K1 = randi(PopSize, 1, N);
    K2 = site2(K1, 2);  % 最近邻索引（第2个）

    index = zeros(1, N);
    for i = 1:N
        if FrontNo(K1(i)) < FrontNo(K2(i))
            index(i) = K1(i);
        elseif FrontNo(K1(i)) == FrontNo(K2(i))
            if SpCrowdDis(K1(i)) >= SpCrowdDis(K2(i))
                index(i) = K1(i);
            else
                index(i) = K2(i);
            end
        else
            index(i) = K2(i);
        end
    end
end

function n = norm_calnum(x,tag)
    if ~exist('x','var') || isempty(x) || ~isnumeric(x)
        fprintf('!!! DEBUG[%s]: calnum is empty or non-numeric: %s\n', tag, mat2str(x));
        n = 0;
    else
        n = double(sum(x(:)));
        if ~isfinite(n) || isnan(n)
            fprintf('!!! DEBUG[%s]: calnum is invalid value (NaN/Inf): %s\n', tag, mat2str(x));
            n = 0;
        end
    end
end