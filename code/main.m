% Please kindly cite the paper Junyi Guan, Sheng li, Jinhui Zhu, Xiongxiong He, and Jiajia Chen 
% "Fast main density peak clustering within relevant regions via a robust decision graph" 
% Pattern Recognition,2024
% The code was written by Junyi Guan in 2023.

clear all;close all;clc;

%% load data  
load data/Agg
data_with_lable = Agg; 

%% deduplicate data
data_x = unique(data_with_lable,'rows');
if size(data_x,1) ~= size(data_with_lable,1)
    data_with_lable = data_x;
end
answer = data_with_lable(:,end);
data = data_with_lable(:,1:end-1);
n = size(data,1);

%% parameters 
k = round(sqrt(n))+1; %default
cv = 1; %default

%% R-MDPC clustering
[CL,rho,delta,centers,runtime] =  R_MDPC(data,cv,k);

%% real centers
re_centers = realcenter(answer,rho);

%% show result
close all
resultshow(data,CL,centers,rho,delta,centers,min(rho),re_centers);

%% evaluation
DGCI = GetDGCI(rho,delta,re_centers,n); %% graph clarity
[AMI,ARI,FMI] = Evaluation(CL,answer); %% clustering result

%% clustering result
result = struct;
result.AMI = AMI;
result.ARI = ARI;
result.FMI = FMI;
result.DGCI = DGCI;
result.k = k;
result.cv = cv;
result.NCLUST = length(centers);
result.runtime = runtime;
result




