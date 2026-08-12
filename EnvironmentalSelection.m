function [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,N)
% The environmental selection of SparseEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% 在 EnvironmentalSelection 函数内部
%disp('--- Debugging Population Variable ---');
%whos Population
%if iscell(Population)
%    disp('Population is a cell array.');
%   if isempty(Population)
%        disp('Population is empty.');
%    else
%        disp(['Size of Population: ', num2str(size(Population))]);
%        disp(['Type of first element: ', class(Population{1})]);
%        % 打印第一个元素的类型，如果是 struct，则说明没问题
%        if isstruct(Population{1})
%            disp('First element is a struct. All good.');
%        else
%            disp('First element is NOT a struct. This is the issue!');
%        end
%    end
%else
%    disp('Population is NOT a cell array. This is the issue!');
%    disp(['Type of Population: ', class(Population)]);
%end
%disp('-------------------------------------');
    %% Delete duplicated solutions
    PopObjs        = extract(Population);
    PopDecs        = extract_d(Population);
    [~,uni]        = unique(PopDecs,'rows');
    if length(uni) == 1
        [~,uni] = unique(PopObjs,'rows');
    end
    Population = Population(uni);
    PopObjs    = PopObjs(uni,:);
    PopDecs    = PopDecs(uni,:);
    Dec        = Dec(uni,:);
    Mask       = Mask(uni,:);

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(PopObjs,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Mask,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    % === 修复：确保请求数量不超过实际数量 ===
    num_to_select = N - sum(Next);
    num_available = length(Last);
    num_to_select = min(num_to_select, num_available);
    Next(Last(Rank(1:num_to_select))) = true;
    %Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
    Dec        = Dec(Next,:);
    Mask       = Mask(Next,:);
end
