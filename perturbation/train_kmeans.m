function centroid = train_kmeans(features,labels,nBoot,ratio)
% train_kmeans unsupervised clustering using the k-means algorithm. Robust
% centroids are obtained by bootstrapping the data such that the proportion 
% of seizure free and non-seizure free patients are equal to the variable
% 'ratio'.
%
%   centroid = train_kmeans(features,labels,nBoot,ratio) applies the
%   k-means algorithm (k=2) to classify the features (Px3 double vector) by
%   randomly sampling total cohort nBoot times such that each sub-sample
%   contains ratio (a number from 0 to 1) seizure-free patients and 
%   non-seizure-free patients. The function also takes in 'labels', which 
%   is a Px1 array, containing information on surgical outcome
%   (seizure-free=1, non-seizure-free=2). The function returns the final 
%   centroid (2x3), obtained by averaging over all the bootstrapped centroids 
%   channel-level features (NxM) and their corresponding coordinates (Nx3),
%   where N is the number of channels, and M is the number of
%   observations. 
%
%   See also computeSPMap and computeCorrelation.

% Selecting the seizure-free and non-seizure-free patients
X_sf = features(labels==1,:);
X_nsf = features(labels==2,:);

engelIACentroid_bs = [];
engelIICentroid_bs = [];

% bootstrap sampling
for i = 1:nBoot
    % Randomly suffling seizure-free and non-seizure-free patients
    idx_sf = randperm(size(X_sf,1))';
    idx_nsf = randperm(size(X_nsf,1))';

    % Obtaining number of patients to sub-sample
    M1 = round(ratio*size(idx_sf,1));
    M2 = round(ratio*size(idx_nsf,1));
    
    % Selecting sub-samples 
    tmp_sf = X_sf(idx_sf,:); train_sf = tmp_sf(1:M1,:); 
    tmp_nsf = X_nsf(idx_nsf,:); train_nsf = tmp_nsf(1:M2,:);

    
    X_train = [train_sf; train_nsf];
    idx_rand = randperm(size(X_train,1));

    % Estimating centroids using k-means with k=2
    [~,centroid] = kmeans(X_train(idx_rand,:),2,'Replicates',5);
    
    % To order the centroids, they are sorted by distance to origin
    dists = pdist2(centroid,[0 0 0]);
    [~,cIdxs] = sort(dists,'ascend'); % good to bad sampling

    % re-ordering centroid accordingly from good to bad
    centroid = centroid(cIdxs,:);

    engelIACentroid_bs(i,:) = centroid(1,:);
    engelIICentroid_bs(i,:) = centroid(2,:);
end

% Mean is used for analysis
engelIACentroid = mean(engelIACentroid_bs);
engelIICentroid = mean(engelIICentroid_bs);

centroid = [engelIACentroid; engelIICentroid];
end