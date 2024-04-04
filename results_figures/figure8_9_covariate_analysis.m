% This script reproduces Figures 8 and 9 of the main article.
%% 
load('SPMap_features.mat'); % Loading SP map features extracted from MNI and CHUGA
addpath('functions');
addpath('demo_data')
addpath('perturbation')
addpath(genpath('results_figures\sp_results'))

% Specify center (1=MNI,2=CHUGA)
center = 1;

% Obtaining MNI covariates
centers = features(:,end);
label_mni = features(centers == 1,5);
X_mni = features(centers == 1,[1 2 4]); %Q1,Q2,Q4
% Column 1: SOZ Volume; Column 2: Resection Volume, Column 3: Resected SOZ Volume, Column 4: Percent SOZ removed; Column 5: Palliative;
cov_mni = covariates(centers == 1,:)/1000; % Converting from mm^3 to cm^3

% Obtaining CHUGA covariates
label_chuga = features(centers == 2,5);
X_chuga = features(centers == 2,[1 2 4]);
cov_chuga = covariates(centers == 2,:)/1000; % Converting from mm^3 to cm^3
%% Setting up centroids as estimated in cluster_analysis.m
centroid = [0.3319    0.1137    0.0696;
            0.4401    0.3843    0.2468];
[~,idx_test_mni] = pdist2(centroid,X_mni,'euclidean','smallest',1); idx_test_mni = idx_test_mni';  % MNI cluster indices
[~,idx_test_chuga] = pdist2(centroid,X_chuga,'euclidean','smallest',1); idx_test_chuga = idx_test_chuga'; % CHUGA cluster indices 
%% Select cut-offs based on clinical covariate (as done in the manuscript) for additional analysis
cutoff_mni = cov_mni(:,3) > 3305/1000;
cutoff_chuga = cov_chuga(:,3) > 3305/1000;
%% Figure 8: Plotting SOZ volumes across surgical outcome and cluster indices
figure;
red =  [255 94 105]/255;
green = [157, 216, 102]/255;
colors = [red;green];

% Figure 8a: Comparing across clusters
Z = cov_mni(:,1);
g = idx_test_mni;  
[p,d,gr]=boxplot_gramm(g,Z,colors,1,1);
gr(1,1).set_names('x','','y','SOZ Volume (cm^3)');
results = [{'Figure 8a'} p d];

% Figure 8b: Comparing within the first cluster across surgical outcome
clusterID = 1; 
Z = cov_mni(idx_test_mni==clusterID,1); 
g = label_mni(idx_test_mni==clusterID); 
[p,d,gr]=boxplot_gramm(g,Z,colors,1,2,gr);
gr(1,2).set_names('x','','y','SOZ Volume (cm^3)');
results(2,:) = [{'Figure 8b'} p d];

% Figure 8c: Correcting for insufficient resections
Z = cov_mni(cutoff_mni,1);
g = idx_test_mni(cutoff_mni); 
[p,d,gr]=boxplot_gramm(g,Z,colors,2,1,gr); 
gr(2,1).set_names('x','','y','SOZ Volume (cm^3)');
results(3,:) = [{'Figure 8c'} p d];

% Figure 8d: Correcting for insufficient resections
Z = cov_chuga(cutoff_chuga,1);
g = idx_test_chuga(cutoff_chuga); 
[p,d,gr]=boxplot_gramm(g,Z,colors,2,2,gr); 
gr(2,2).set_names('x','','y','SOZ Volume (cm^3)');
results(4,:) = [{'Figure 8d'} p d];

gr.draw();

set(gcf,"Position",[ 21   189   905   747]);
movegui('west');
%% Figure 9: Correlating percent SOZ removed with distance to c1
figure;
colors = [green;red];
% Computing distances from each patient in the feature space to c1
[distances_c1] = pdist2(centroid(1,:),X_mni,'euclidean')';

% Figure 9a: Correlation between proportion of SOZ removed and distance
% from c1 (seizure-free)
Z = cov_mni(label_mni==1,4); % Percent SOZ resected
g = distances_c1(label_mni==1);
[pval,rho,gr]=scatter_gramm(g,Z,idx_test_mni(label_mni==1),colors,green,2,1,1);
gr(1,1).set_names('x','Distance from c_1','y','SOZ removed (%)');
results(5,:) = [{'Figure 9a'} pval rho];

% Figure 9b: Correlation between proportion of SOZ removed and distance
% from c1 (non-seizure-free)
Z = cov_mni(label_mni==2,4); % Percent SOZ resected
g = distances_c1(label_mni==2);
[pval,rho,gr]=scatter_gramm(g,Z,idx_test_mni(label_mni==2),colors,red,1,1,2,gr);
gr(1,2).set_names('x','Distance from c_1','y','SOZ removed (%)');
results(6,:) = [{'Figure 9b'} pval rho];

% Figure 9c: Correlation between resection volume and distance
% from c1 (seizure-free)
Z = cov_mni(label_mni==1,2); % Resection volume
g = distances_c1(label_mni==1);
[pval,rho,gr]=scatter_gramm(g,Z,idx_test_mni(label_mni==1),colors,green,2,2,1,gr);
gr(2,1).set_names('x','Distance from c_1','y','Resected Volume (cm^3)');
results(7,:) = [{'Figure 9c'} pval rho];

% Figure 9d: Correlation between resection volume and distance
% from c1 (seizure-free)
Z = cov_mni(label_mni==2,2); % Resection volume
g = distances_c1(label_mni==2);
[pval,rho,gr]=scatter_gramm(g,Z,idx_test_mni(label_mni==2),colors,red,1,2,2,gr);
gr(2,2).set_names('x','Distance from c_1','y','Resected Volume (cm^3)');
results(8,:) = [{'Figure 9d'} pval rho];

gr.draw();
set(gcf,'Position',[975   133   938   802]);
movegui('east')
