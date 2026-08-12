function PopObj        = extract(Population)
%% PopObj=Population.objs
        PopObj=zeros(size(Population,2),2);
        for i=1:size(PopObj,1)
            PopObj(i,:)=Population{1,i}.obj;
        end
end