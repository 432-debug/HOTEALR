function [CSN] = spcc_method(data,i)
%construct sample specific network

mean_a=mean(data,2);
std_a=var(data,0,2);
    v=data(:,i);
    x=(v-mean_a)./std_a;
    x(isnan(x))=0;
    x(abs(x)==inf)=0;
    
    cand1=(x*x');
    x=cand1;
    xx=triu(x);
    a_x=abs(xx(xx~=0));
   
    CSN=xx;




end

