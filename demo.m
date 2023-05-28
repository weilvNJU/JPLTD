clear; close all;
tic
addpath E:\聚类论文实验2\代码整理\ClusteringMeasure
addpath E:\聚类论文实验2\代码整理\twist
path = 'E:\聚类论文实验2\代码整理\data\Index\';

%%
load('E:\聚类论文实验2\代码整理\data\ORL.mat');
gt = Y;
name = 'ORL';
percentDel = 0.1;
Datafold= [path,'Index_',name,'_percentDel_',num2str(percentDel),'.mat'];
load(Datafold)
param.lambda = 5;
param.theta = 0.1;
k = 100;
cls_num = numel(unique(Y));
%%
fPCA = 0;
perf = [];
for idx = 1
    Xc = X;
    ind = Index{idx};
    for i=1:length(Xc)
        Xci = Xc{i};
        indi = ind(:,i);
        pos = find(indi==0);
        Xci(:,pos)=[];
        Xc{i} = Xci;
    end
    [G, B, P, Loss] = JPLTD(Xc, ind, Y, param, k);
    for rp = 1:20
        [Clus] = SpectralClustering(G, cls_num);
        result =  Clustering8Measure(gt,Clus);
        perf  = [perf ; result*100];
    end
end
mean(perf)
std(perf)
perf = [];
toc


