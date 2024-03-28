
% Please kindly cite the paper Junyi Guan, Sheng li, Jinhui Zhu, Xiongxiong He, and Jiajia Chen
% "Fast main density peak clustering within relevant regions via a robust decision graph"
% Pattern Recognition,2024
% The code was written by Junyi Guan in 2023.

function [CL,rho,delta,centers,runtime] = R_MDPC(data,cv,k)
fprintf('R-MDPC clusteing! :)\n')

import java.util.LinkedList
import Library.*
close all
tic;

%% Normalization
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;

[~,d]  = size(data); % d:dimensions
k_b = round(k/2)+1;

%% Fast KNN calculation based on Kd-tree (when dimension is not large than 10)
if d<=10
    [knn_idx,knn_dist] = knnsearch(data,data,'k',2*k);
else
    dist = squareform(pdist(data));
    [knn_dist,knn_idx] = sort(dist,2);
end

%% Obtain a suitable density distribution
lambda = 0;
for i = 1:100
    lambda = 0.1+lambda;
    rho=sum(knn_dist(:,2:k).^lambda,2).^-1; %% KNN(+)-Density
    current_cv = std(rho)/mean(rho); %% ±äÒìÏµÊý
    if  current_cv> cv
        break
    end
end

time0 = toc;

tic;


%% single-peak clustering
[peaks,rho_peaks,sub_label,n_peaks,ordrho,parent] = singlePeakClustering(rho,knn_idx,k);

%% region construction
[region_peaks,region_label,n_regions,jun] = regionConstruction(rho,sub_label,peaks,knn_idx,k_b);

%% delta and rho calculation within regions
[delta,rho,parent]= deltaRhoCalculationWithihRegions(data,region_label,region_peaks,rho,peaks,parent,knn_idx,knn_dist);

%% satellite peak attenuation
[delta] = attenuartor(delta,jun,peaks,rho_peaks);
time1 = toc;

%% decisionAllocation
[CL,NCLUST,centers,time2] =  decisionAllocation(rho,delta,peaks,n_regions,ordrho,parent);
runtime = time0 + time1 + time2;



%% functions

%%% single-peak clustering function
function[peaks,rho_peaks,sub_label,n_peaks,ordrho,parent] = singlePeakClustering(rho,knn_idx,k)
n = length(rho);
[~,ordrho]=sort(rho,'descend');
peaks = []; %peaks: density peaks
parent = zeros(1,n);%parent node
for i=1:n
    depth(ordrho(i))= 0;
    for j=2:k
        neigh=knn_idx(ordrho(i),j);
        if(rho(ordrho(i))<rho(neigh))
            parent(ordrho(i))=neigh;
            depth(ordrho(i))= depth(neigh) + 1;
            break
        end
    end
    if depth(ordrho(i))==0
        peaks = [ordrho(i) peaks];
    end
end
%% label initialization
for i=1:n
    sub_label(i)=-1;
end

n_peaks = length(peaks);%% n_peaks:Number of root nodes
sub_label(peaks) = (1:n_peaks);
rho_peaks = rho(peaks);
for i=1:n
    if (sub_label(ordrho(i))==-1)
        sub_label(ordrho(i))=sub_label(parent(ordrho(i)));
    end
end

%%% region construction function
function [region_peaks,region_label,n_regions,jun] = regionConstruction(rho,sub_label,peaks,knn_idx,k1)
n_peaks = length(peaks);
n = length(sub_label);
connect_matrix = zeros(n_peaks,n_peaks);% connect information between sub-clusters
jun = zeros(n_peaks,1);

for i=1:n
    for j = 2:k1
        i_neigr = knn_idx(i,j);
        if (sub_label(i)~=sub_label(i_neigr)) && (ismember(i,knn_idx(i_neigr,2:k1)))
            if rho(i_neigr)>jun(sub_label(i)) && rho(peaks(sub_label(i)))< rho(peaks(sub_label(i_neigr)))
                jun(sub_label(i)) = rho(i_neigr);
            end
            connect_matrix(sub_label(i),sub_label(i_neigr)) = 1;
        end
    end
end

zeros(n_peaks,n_peaks);
sim_list = [];
for i = 1:n_peaks-1
    for j = i+1:n_peaks
        if connect_matrix(i,j) && connect_matrix(j,i)
            sim = 1;
        else
            sim = 0;
        end
        sim_list = [sim_list sim];
    end
end
%% SingleLink clustering of sub-clusters according to SIM
SingleLink = linkage(1-sim_list,'single');
sim_inf = 1.-SingleLink(:,3);
n_regions =  n_peaks-length(find(sim_inf>0));

if n_regions>1
    region_label_peak = cluster(SingleLink,n_regions);
else
    region_label_peak = ones(n,1);
end
%% assign region label
for i=1:n_peaks
    region_label(sub_label==i) = region_label_peak(i); %% region_label: the region label
end
%% find region peaks
for i=1:n_regions
    peaks_within_region = peaks(region_label(peaks)==i);
    rhomax_peak_index = find(rho(peaks_within_region)==max(rho(peaks_within_region)));
    region_peaks(i) = peaks_within_region(rhomax_peak_index(1)) ;
end

%%%
function [delta,rho,parent]= deltaRhoCalculationWithihRegions(data,region_label,region_peaks,rho,peaks,parent,knn_idx,knn_dist)
rho_peaks = rho(peaks);
n_peaks = length(peaks);
n_regions = length(region_peaks);
[n,d] = size(data);
delta = zeros(n,1);%% delta initialization
delta(peaks)=Inf;%% delta initialization

%% the distance matrix between RNs and its region points
tic;


if d<=10
    min_peakRho = min(rho_peaks);
    largeRho_idx = find(rho>=min_peakRho);
    n_larg = length(largeRho_idx);
    data_larg = data(largeRho_idx,:);
    [idx_peaks,dist_peaks] = knnsearch(data_larg,data(peaks,:),'k',n_larg); %%
else
    largeRho_idx = (1:n);
    idx_peaks  = knn_idx(peaks,:);
    dist_peaks = knn_dist(peaks,:);
end

for i=1:n_peaks
    p1 = peaks(i);
    idx_peaks_1_all = idx_peaks(i,:);
    if d<=10
        idx_peaks_1_R = idx_peaks(i, region_label(largeRho_idx(idx_peaks_1_all))==region_label(p1));
        dist_peaks_1_R = dist_peaks(i,region_label(largeRho_idx(idx_peaks_1_all))==region_label(p1));
    else
        idx_peaks_1_R = idx_peaks(i, region_label(largeRho_idx(idx_peaks_1_all))==region_label(p1));
        dist_peaks_1_R = dist_peaks(i,region_label(largeRho_idx(idx_peaks_1_all))==region_label(p1));
    end
    
    if ~ismember(p1,region_peaks)
        for j = 1:length(idx_peaks_1_R)
            PN_p1=largeRho_idx(idx_peaks_1_R(j));
            if(rho(p1)<rho(PN_p1))
                delta(p1)=dist_peaks_1_R(j);
                parent(p1)=PN_p1;
                break
            end
        end
    else
        delta(p1) = 0;
    end
    
end
delta(region_peaks) = max(delta)*1.5;

if n_regions==1
    delta(region_peaks)=0;
    delta(region_peaks) = max(delta)*1.5;
end

for i=1:n_regions
    rho(region_label==i) = rho(region_label==i)./rho(region_peaks(i));
end

%%% attenuartor function
function [delta] = attenuartor(delta,jun,peaks,rho_peaks)
n_peaks = length(peaks);
for i=1:n_peaks
    phi = min(1,jun(i)/rho_peaks(i));
    attenuation = (1-phi);
    delta(peaks(i)) =   delta(peaks(i)) * attenuation;
end

%%% decision and allocation
function [CL,NCLUST,centers,time2] =  decisionAllocation(rho,delta,peaks,n_regions,ordrho,parent)
n = length(rho);
n_peaks = length(peaks);

figure;
plot(rho,delta,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
%     plot(rho(r_centers),delta(r_centers),'o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
text(0.9,1.2,[num2str(n_regions),' regions'],'FontSize',15);
hold on
title ('Decision Graph','FontSize',15.0);
xlabel ('\rho');
ylabel ('\delta');
axis([min(rho) max(rho) min(delta) max(delta)]);
rect = getrect;
rhomin=rect(1);
deltamin=rect(2);


%% allocation
tic
for i=1:n
    CL(i)=-1;
end
NCLUST=0;


for i=1:n_peaks
    if rho(peaks(i))>rhomin && delta(peaks(i))>deltamin
        NCLUST=NCLUST+1;
        CL(peaks(i))=NCLUST;
        centers(NCLUST)=peaks(i);
    end
end


for i=1:n
    if (CL(ordrho(i))==-1)
        CL(ordrho(i))=CL(parent(ordrho(i)));
    end
end

time2 = toc;






