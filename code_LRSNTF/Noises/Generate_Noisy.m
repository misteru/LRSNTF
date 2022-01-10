addpath(genpath('Data'));
addpath(genpath('Noises'));
addpath(genpath('lib'));

load('Ori_WDC.mat');
load('WDC_level.mat');
load('WDC_band.mat');
%%
Ohsi = Img;

if max(Ohsi(:))>1
    Ohsi=my_normalized(Ohsi);
end

Nway = size(Ohsi);

%% noise_case1: G:0.1-0.2
Nhsi = zeros(Nway);
for j=1:Nway(3)
    Nhsi(:,:,j) = Ohsi(:,:,j)+sigma_n3(j)*randn(Nway(1),Nway(2));
end
save('./Noises/WDC_case1','Ohsi','Nhsi'); % save results  


%% noise_case2 G:0.1-0.2 S;0.1-0.2
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end
save('./Noises/WDC_case2','Ohsi','Nhsi'); % save results   


%% noise_case3 G:0.1-0.2 S;0.1-0.2 stripes: 30%的band, 每个band为6-15条 
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end

for i = 1:length(add_band)
    stripenum = randperm(10,1)+5;
    locolumn    = randperm(Nway(2),stripenum);
    Nhsi(:,locolumn,add_band(i))=0.2*rand(1)+0.6;
end
save('./Noises/WDC_case3','Ohsi','Nhsi'); % save results  



%% noise_case4 G:0.1-0.2 S;0.1-0.2 deadline: 30%的band, 每个band为6-10条， 宽度为1-3
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end

for i = 1:length(add_band)
    deadlinenum = randperm(5,1)+5;
    locolumn    = randperm(Nway(2)-2,deadlinenum);
    an          = funrand(3,deadlinenum);
    loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
    Nhsi(:,locolumn(loc1),add_band(i))=0; 
    Nhsi(:,locolumn(loc2),add_band(i))=0;
    Nhsi(:,locolumn(loc2)+1,add_band(i))=0;
    Nhsi(:,locolumn(loc3),add_band(i))=0;
    Nhsi(:,locolumn(loc3)+1,add_band(i))=0; 
    Nhsi(:,locolumn(loc3)+2,add_band(i))=0; 
end
save('./Noises/WDC_case4','Ohsi','Nhsi'); % save results  



%% noise_case5 G:0.1-0.2 S;0.1-0.2  stripes: 30%的band, 每个band为6-15条  deadline: 30%的band, 每个band为5-10条， 宽度为1-3
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end

for i = 1:length(add_band)
    stripenum = randperm(10,1)+5;
    locolumn    = randperm(Nway(2),stripenum);
    Nhsi(:,locolumn,add_band(i))=0.2*rand(1)+0.6;
end

for i = 1:length(add_band)
    deadlinenum = randperm(5,1)+5;
    locolumn    = randperm(Nway(2)-2,deadlinenum);
    an          = funrand(3,deadlinenum);
    loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
    Nhsi(:,locolumn(loc1),add_band(i))=0; 
    Nhsi(:,locolumn(loc2),add_band(i))=0;
    Nhsi(:,locolumn(loc2)+1,add_band(i))=0;
    Nhsi(:,locolumn(loc3),add_band(i))=0;
    Nhsi(:,locolumn(loc3)+1,add_band(i))=0; 
    Nhsi(:,locolumn(loc3)+2,add_band(i))=0; 
end
save('./Noises/WDC_case5','Ohsi','Nhsi'); % save results   