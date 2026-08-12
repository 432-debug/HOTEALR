function [rbm,allOne,other,Site] = trainRBM(Mask,UpdateRatio,N,current_subspace)

    allZero = all(~Mask,1);
    allOne  = all(Mask,1);
    other   = ~allZero & ~allOne;
    other_idx = find(other);  % 将other逻辑数组转为索引数组（哪些列是目标列）
    other_sub = intersect(other_idx, current_subspace);% other_sub：current_subspace中属于other的列（索引数组，列数=目标维度数）

    Site1 = false(1,floor(N/2));
    Site2 = false(1,floor(N/2));
    Site = [Site1,Site2];
    
    rbm = [];
    
    if rand < UpdateRatio %随机数小于后代存活率才使用RBM
        
        Site1(randperm(floor(N/2),randi(ceil(N/2*0.1)))) = true; % 从N中随机选择10%个体使用RBM生成子代
        
        Site2 = Site1;

        Site = [Site1,Site2];

        %K = ceil(sum(mean(abs(Mask(:,other).*Dec(:,other))>1e-6,1)));
        K = ceil(sum(mean(abs(Mask(:,other_sub))>1e-6,1)));
        
        K = min(K,ceil(size(Mask,2)*0.2));
        
        %rbm = RBM(sum(other),K,8,1,0,0.4,0.15);
        rbm = RBM(length(other_sub),K,10,1,0,0.5,0.1);
        
        rbm.train(Mask(:,other_sub));

    end
end



% function [rbm,allOne,other,Site] = trainRBM(Mask,UpdateRatio,N,current_subspace)
% 
%     allZero = all(~Mask,1);
%     allOne  = all(Mask,1);
%     other   = ~allZero & ~allOne;
%     other_idx = find(other);  % 将other逻辑数组转为索引数组（哪些列是目标列）
%     % other_sub = intersect(other_idx, current_subspace);% other_sub：current_subspace中属于other的列（索引数组，列数=目标维度数）
% 
%     Site1 = false(1,floor(N/2));
%     Site2 = false(1,floor(N/2));
%     Site = [Site1,Site2];
% 
%     rbm = [];
% 
%     if rand < UpdateRatio %随机数小于后代存活率才使用RBM
% 
%         Site1(randperm(floor(N/2),randi(ceil(N/2*0.1)))) = true; % 从N中随机选择10%个体使用RBM生成子代
% 
%         Site2 = Site1;
% 
%         Site = [Site1,Site2];
% 
%         %K = ceil(sum(mean(abs(Mask(:,other).*Dec(:,other))>1e-6,1)));
%         K = ceil(sum(mean(abs(Mask(:,other))>1e-6,1)));
% 
%         K = min(K,ceil(size(Mask,2)*0.2));
% 
%         %rbm = RBM(sum(other),K,8,1,0,0.4,0.15);
%         rbm = RBM(sum(other),K,10,1,0,0.5,0.1);
% 
%         rbm.train(Mask(:,other));
% 
%     end
% end

