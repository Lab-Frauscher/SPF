% Loading SP map features to find the optimal selection 
addpath(genpath('results_figures'))
load('SPMap_features.mat');

% Training is on the MNI cohort only
centers = features(:,end);
Z = features(centers == 1,1:4);

%% Feature selection (Supplementary Table S2)
% Computing cluster compactness and seperabiity for all combinations
rng('default'); % for reproducability
dunns_index = [];
combinations = {[]};
for k = 1:9 % For each numer of clusters
    index = 1;
    for i = 1:4
        for j = i+1:4 
            tmp = Z; tmp(:,[i j]) = [];
            indices = setdiff(1:4,[i j]);
            combinations(index,1) = {"Q"+ string(indices(1)) + ",Q" + string(indices(2))};
            distM=squareform(pdist(tmp));
            [ids]=kmeans(tmp,k+1);
            dunns_index(index,k) = dunns(k,distM,ids);
            index = index + 1;
        end
    end
end
for k = 1:9
    index2 = index;
    for i = 1:4
        tmp = Z; tmp(:,i) = [];
        indices = setdiff(1:4,i);
        combinations(index2,1) = {"Q"+ string(indices(1)) + ",Q" + string(indices(2)) + ",Q" + string(indices(3))};
        distM=squareform(pdist(tmp));
        [ids]=kmeans(tmp,k+1);
        dunns_index(index2,k) = dunns(k,distM,ids);
        index2 = index2 + 1;
    end
end
table = cat(2,combinations,num2cell(dunns_index));
disp(table)
[~,idx]=max(dunns_index,[],1); % Obtaining optimal combination for each k
disp('Optimal feature combinations for each k')
disp(strjoin("〈" + string(combinations(idx))' + "〉",''))
%% Optimal number of clusters
f_opt = Z(:,[1 2 4]); % Most frequenct combination with high dunn's index
ks = 1:10;
error = [];
for k=ks
    [idx_test,centroid,sumd] = kmeans(f_opt,k,'Replicates',5);
    error = [error; sum(sumd)];
end

ii = knee_pt(error,ks);
plot(ks,error,'k-','LineWidth',1); hold on;
plot(ks(ii),error(ii),'rx','LineWidth',2,'MarkerSize',15);
xlabel('Number of clusters');
ylabel('e_d(m)','Interpreter','tex')
set(gca,'FontSize',20)
xlim([ks(1) ks(end)])
ylim([min(error) max(error)])