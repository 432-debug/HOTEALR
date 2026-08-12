function PopDec = extract_d(Population)
    % 增加判空逻辑
    if isempty(Population)
        PopDec = []; 
        return;
    end
    
    % 获取种群大小
    N = length(Population);
    % 安全获取第一个个体的维度
    if N > 0 && isfield(Population{1}, 'dec')
        D = size(Population{1}.dec, 2);
        PopDec = zeros(N, D);
        for i = 1:N
            PopDec(i, :) = Population{i}.dec;
        end
    else
        PopDec = [];
    end
end
%function PopObj        = extract_d(Population)
%% PopObj=Population.decs
%        PopObj=zeros(size(Population,2),size(Population{1, 1}.dec,2));
%        for i=1:size(PopObj,1)
%            PopObj(i,:)=Population{1,i}.dec;
%        end
%end