clc
clear
a1=importdata('mouse_scRNAseq_CD45negative_cells_0W_NOR.csv');
a2=importdata('mouse_scRNAseq_CD45negative_cells_12W_IFN.csv');
a3=importdata('mouse_scRNAseq_CD45negative_cells_22W_DYS.csv');
a4=importdata('mouse_scRNAseq_CD45negative_cells_24W_CIS.csv');
a5=importdata('mouse_scRNAseq_CD45negative_cells_26W_ICA.csv');

%1-5581
%5582-12221
%12222-18483
%18484-22422
%22423-27691

data1=a1.data;
data2=a2.data;
data3=a3.data;
data4=a4.data;
data5=a5.data;

a=importdata('mouse.xlsx');
Net=upper(a.textdata);

genes=upper(a1.textdata(2:end,1));

 
Final_genes=intersect(unique(Net),genes);%基因名字


[z1,z2]=ismember(Net,Final_genes);
z=z1(:,1).*z2(:,2);
New_Net=z2(find(z~=0),:);

[numID,~]=ismember(Final_genes,genes);

A=zeros(length(Final_genes));
for i=1:size(New_Net,1)
    
    A(New_Net(i,1),New_Net(i,2))=1;
    A(New_Net(i,2),New_Net(i,1))=1;
     
end

%*******************Data1*******************************
New_data1=data1(numID,:);%参数1
sub_network1=[];

for j=1:size(New_data1,2)
    
    j
    tic
    csn = csnet(New_data1,j);
    spccsn = spcc_method(New_data1,j);
    Sample_net=csn{1,j}.*spccsn.*A;
    Sample_net(isnan(Sample_net))=0;
    [sub_z1,sub_z2]=find(triu(Sample_net~=0));
    sub_score=[];
    for i=1:length(sub_z1)
        sub_score(i,1)=Sample_net(sub_z1(i,1),sub_z2(i,1));
    end
    sub_network1{j,1}=[sub_z1 sub_z2 sub_score];%结果
    toc

end

%*******************Data2*******************************
New_data2=data2(numID,:);%参数1
sub_network2=[];

for j=1:size(New_data2,2)
    
    j
    tic
    csn = csnet(New_data2,j);
    spccsn = spcc_method(New_data2,j);
    Sample_net=csn{1,j}.*spccsn.*A;
    Sample_net(isnan(Sample_net))=0;
    [sub_z1,sub_z2]=find(triu(Sample_net~=0));
    sub_score=[];
    for i=1:length(sub_z1)
        sub_score(i,1)=Sample_net(sub_z1(i,1),sub_z2(i,1));
    end
    sub_network2{j,1}=[sub_z1 sub_z2 sub_score];%结果
    toc

end

%*******************Data3*******************************
New_data3=data3(numID,:);%参数1
sub_network3=[];

for j=1:size(New_data3,2)
    
    j
    tic
    csn = csnet(New_data3,j);
    spccsn = spcc_method(New_data3,j);
    Sample_net=csn{1,j}.*spccsn.*A;
    Sample_net(isnan(Sample_net))=0;
    [sub_z1,sub_z2]=find(triu(Sample_net~=0));
    sub_score=[];
    for i=1:length(sub_z1)
        sub_score(i,1)=Sample_net(sub_z1(i,1),sub_z2(i,1));
    end
    sub_network3{j,1}=[sub_z1 sub_z2 sub_score];%结果
    toc

end


%*******************Data4*******************************
New_data4=data4(numID,:);%参数1
sub_network4=[];

for j=1:size(New_data4,2)
    
    j
    tic
    csn = csnet(New_data4,j);
    spccsn = spcc_method(New_data3,j);
    Sample_net=csn{1,j}.*spccsn.*A;
    Sample_net(isnan(Sample_net))=0;
    [sub_z1,sub_z2]=find(triu(Sample_net~=0));
    sub_score=[];
    for i=1:length(sub_z1)
        sub_score(i,1)=Sample_net(sub_z1(i,1),sub_z2(i,1));
    end
    sub_network4{j,1}=[sub_z1 sub_z2 sub_score];%结果
    toc

end


%*******************Data5*******************************
New_data5=data5(numID,:);%参数1
sub_network5=[];

for j=1:size(New_data5,2)
    
    j
    tic
    csn = csnet(New_data5,j);
    spccsn = spcc_method(New_data5,j);
    Sample_net=csn{1,j}.*spccsn.*A;
    Sample_net(isnan(Sample_net))=0;
    [sub_z1,sub_z2]=find(triu(Sample_net~=0));
    sub_score=[];
    for i=1:length(sub_z1)
        sub_score(i,1)=Sample_net(sub_z1(i,1),sub_z2(i,1));
    end
    sub_network5{j,1}=[sub_z1 sub_z2 sub_score];%结果
    toc

end


%[P1,R1]=corrcoef(data1');
%R1(R1>0.05)=0;
%N1=P1.*R1.*A;

%[P2,R1]=corrcoef(data2');
%R1(R1>0.05)=0;
%N2=P2.*R1.*A;


%[P3,R1]=corrcoef(data3');
%R1(R1>0.05)=0;
%N3=P3.*R1.*A;

%[P4,R1]=corrcoef(data4');
%R1(R1>0.05)=0;
%N4=P4.*R1.*A;


%[P5,R1]=corrcoef(data5');
%R1(R1>0.05)=0;
%N5=P5.*R1.*A;

 
   name='mouse_scRNAseq_data';
   samplename=strcat(name,'.mat');
   save(samplename,'sub_network1','sub_network2','sub_network3','sub_network4','sub_network5','Final_genes')

    




