function [Population,Dec,Mask,FrontNo,SpCrowdDis,sRatio] = EnvironmentalSelection_adr(Population,Dec,Mask,N,len,num)
% The environmental selection of MOEA/PSL

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Delete duplicated solutions
success = false(1,length(Population));
PopObjs        = extract(Population);
PopDecs        = extract_d(Population);
[~,uni] = unique(PopDecs,'rows');
if length(uni) == 1
    [~,uni] = unique(PopObjs,'rows');
end
Population = Population(uni);
PopObjs    = PopObjs(uni,:);
PopDecs    = PopDecs(uni,:);
Dec        = Dec(uni,:);
Mask       = Mask(uni,:);

N          = min(N,length(Population));

%% Non-dominated sorting
[FrontNo,SpCrowdDis] = M_non_domination_scd_sort(Mask,PopObjs);
fnum=0;                                                                 %当前前沿面
while numel(FrontNo,FrontNo<=fnum+1)<=N                      %numel(A,A<5)判断矩阵A中小于5的个数，判断前多少个面的个体能完全放入外部存档
    if numel(FrontNo,FrontNo<=fnum+1)==numel(FrontNo,FrontNo<=fnum)
        fnum=fnum-1;
        break
    end
    fnum=fnum+1;
end
MaxFNo=fnum+1;
%     MaxFNo=max(FrontNo);
%     [FrontNo,MaxFNo] = NDSort(PopObjs,N);
Next = FrontNo < MaxFNo;

%% Calculate the crowding distance of each solution
%     CrowdDis = CrowdingDistance(PopObjs,FrontNo);

%% Select the solutions in the last front based on their crowding distances
Last     = find(FrontNo==MaxFNo);
[~,Rank] = sort(SpCrowdDis(Last),'descend');
Next(Last(Rank(1:N-sum(Next)))) = true;

%% Calculate the ratio of successful offsprings
success(uni(Next)) = true;
s1     = sum(success(len+1:len+num));
s2     = sum(success(len+num+1:end));
sRatio = (s1+1e-6)./(s1+s2+1e-6);
sRatio = min(max(sRatio,0.1),0.9);

% s1     = sum(success(len+1:min(len+num,size(success,2))));
% s2     = sum(success(len+num+1:end));
% 
% UpdateRatio = (s1+s2+1e-6)/(size(success,2)-len+1e-6);
% 
% if num > 1e-6  && s1 > 1e-6
%    UpdateRatio = 0.5 + s1/(s1+s2);
% end

%% Population for next generation
Population = Population(Next);
FrontNo    = FrontNo(Next);
SpCrowdDis = SpCrowdDis(Next);
Dec        = Dec(Next,:);
Mask       = Mask(Next,:);
end