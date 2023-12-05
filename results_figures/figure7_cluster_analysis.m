% This script reproduces Figure 7 in the results of the manuscript.
% It estimates the centroids and the probability model using the 
% bootstrapping method described in the Method section of the article.
%
% The first part of this script loads the features extracted from the SP
% maps of all patients in our dataset. If the center is specificed as MNI,
% then unsupervised clustering is applied to obtain the centroids.
%
% The next part of the script produces a probability model using distances
% from each patient in the feature space to both centroids. The sample
% standard deviation and means are computed and used to calcalate the
% probability that the patient is in the good sampling cluster, and is not
% in the poor sampling cluster. Results are shown for the classification
% performance of the k-means clustering algorithm, and the probability
% model.
%
% The final part of the script will apply the trained model after applying
% a correction on the resected SOZ volume. Only patients with a 
% sufficiently large resected SOZ volume were considered to have complete 
% resections of the SOZ.
%
% The probability models and scatter plots shown in the manuscripts can be
% reproduced by running the code and specifying the appropriate center
addpath('functions');
addpath('perturbation')
addpath('functions')
addpath(genpath('results_figures\sp_results'))
%% 
load('SPMap_features.mat');

center_names = ["MNI";"CHUGA"];
% Specify center (1=MNI,2=CHUGA)
center = 1; %r Uncomment if running directly this script

% Specify training (=1) or testing (=0)
train = center == 1; % MNI is training

centers = features(:,end);
label = features(centers == center,5);
% Q1,Q2,Q4 were found to be the optimal combination
X = features(centers == center,[1 2 4]); 
% Column 1: SOZ Volume; Column 2: Resection Volume, Column 3: Resected SOZ Volume, Column 4: Percent SOZ removed; Column 5: Palliative];
cov = covariates(centers == center,:); 

%% Obtaining threshold on resected SOZ volume to correct for palliative surgeries
% Clinician marked whether patient had an incomplete resection due to
% eloquent cortex
isComplete = covariates(centers==1,end); % MNI data
resectedSOZ = covariates(centers==1,3);
threshold = prctile(resectedSOZ(~isComplete),75); % Upper quadrant on patients marked as palliative

% Select cut-offs based on clinical covariate (as done in the manuscript)
sozResected = cov(:,3);
cutoff = sozResected > threshold; % Applying threshold to correct for insufficient resections of the SOZ
X_c = X(cutoff,:);
label_c = label(cutoff);
cov_c = cov(cutoff,:);
%% K-means algorithm with bootstrapping
% Paramters used for unsupervised clustering
nBoot = 1000;
ratio = 0.75;

rng('default')
if train == 1
    % Estimating robust cluster centroids by bootstrapping MNI data 
    centroid = train_kmeans(X,label,nBoot,ratio);
else % if training option is not chosen, use the centroid estimated using MNI (center=1)
    centroid = [0.3319    0.1137    0.0696;
                0.4401    0.3843    0.2468];
end
%% Figure 7a,b: Performance metrics for the clustering algorithm
[~,idx_test] = pdist2(centroid,X,'euclidean','Smallest',1); idx_test = idx_test';
confusionMatrix = confusionmat(label,idx_test);
TP = confusionMatrix(1,1);
TN = confusionMatrix(2,2);
FP = confusionMatrix(2,1);
FN = confusionMatrix(1,2);

sensitivity = TP / (TP + FN);
specificity = TN / (TN + FP);

gr = plotFeatureSpace(X,label,centroid); % Plotting feature space
%% Estimating parameters of Student's T distribution for data in clusters 1 and 2
[distances_c1] = pdist2(centroid(1,:),X,'euclidean')';
[distances_c2] = pdist2(centroid(2,:),X,'euclidean')'; 
[~,idx_test] = pdist2(centroid,X,'euclidean','smallest',1); idx_test = idx_test'; 

if train == 1    
    mu_1 = mean(distances_c1(idx_test==1));  
    mu_2 = mean(distances_c2(idx_test==2)); 
    
    sigma_1 = std(distances_c1(idx_test==1),1); 
    sigma_2 = std(distances_c2(idx_test==2),1); 
    v_1 = sum(idx_test==1)-1; 
    v_2 = sum(idx_test==2)-1;
else
    mu_1 =  0.1202;
    mu_2 = 0.1242;
    sigma_1 = 0.0485;
    sigma_2 =  0.0604;
    v_1 = 23;
    v_2 = 25;
end
%% Computing p-values and averaging to obtain probabilities
tdist2T = @(t,v) (1-betainc(v./(v+t.^2),v/2,0.5)); 

t_score_c1 = (distances_c1 - mu_1)/(sigma_1);
t_score_c2 = (distances_c2 - mu_2)/(sigma_2);

probability_sampled_c1 = [];
probability_sampled_c2 = [];
for j = 1:length(distances_c1)
    % Probability using cluster 1
    probability_sampled_c1(j,1) = tdist2T(t_score_c1(j),v_1);
    % Probability using cluster 2
    probability_sampled_c2(j,1) = tdist2T(t_score_c2(j),v_2);
end
probability_sampled = 1/2.*(1-probability_sampled_c1) + 1/2.*probability_sampled_c2;
%% Figure 7: Correlating probability model with surgical outcome (need to switch centers to obtain Figure 7c,d)
red =  [255 94 105]/255;
green = [157, 216, 102]/255;

% Figure 7b/c: Probability model and surgical outcomes
Z = probability_sampled;
g = label;

% figure;
% clear gr;
gr(1,2) = gramm('x',ones(size(Z)),'y',Z,'color',g);
gr(1,2).stat_boxplot('width',0.3,'notch',true);
gr(1,2).geom_jitter('alpha',0.5,'dodge',0.7,'width',0.1);
gr(1,2).set_point_options('base_size',15);
gr(1,2).set_text_options('base_size',20);
gr(1,2).set_color_options('map',[green;red],'n_color',2,'n_lightness',1);
gr(1,2).set_names('x','','y','Probability of SOZ sampled');
gr(1,2).axe_property('YLim',[0 1],'XTick',[1 2],'XTickLabel',[{'Seizure-free'};{'Non-seizure-free'}]);
gr(1,2).no_legend();
set(gcf,'units','normalized','outerposition',[0 0 0.3 0.5]);
movegui('center');


% Figure 7d/e: Probability model and surgical outcomes after correction
Z_c = probability_sampled(cutoff);
g_c = label(cutoff);

gr(1,3) = gramm('x',ones(size(Z_c)),'y',Z_c,'color',g_c);
gr(1,3).stat_boxplot('width',0.3,'notch',true);
gr(1,3).geom_jitter('alpha',0.5,'dodge',0.7,'width',0.1);
gr(1,3).set_point_options('base_size',15);
gr(1,3).set_text_options('base_size',20);
gr(1,3).set_color_options('map',[green;red],'n_color',2,'n_lightness',1);
gr(1,3).set_names('x','','y','Probability of SOZ sampled');
gr(1,3).axe_property('YLim',[0 1],'XTick',[1 2],'XTickLabel',[{'Seizure-free'};{'Non-seizure-free'}]);
gr(1,3).no_legend();
set(gcf,'units','normalized','outerposition',[0 0 0.3 0.5]);
movegui('center');

gr.set_title(center_names(center));
gr.draw();
set(gcf,'Position',[ 209,302,1574,652]);
%% Computing performance metrics
[ROC_X,ROC_Y,~,AUC] = perfcurve(g,Z,1,'TVals',0:0.001:1);
[ROC_X_c,ROC_Y_c,~,AUC_c] = perfcurve(g_c,Z_c,1,'TVals',0:0.001:1);

% Before correction statistics
if kstest(Z(g==1)) || kstest(Z(g==2))
    p = ranksum(Z(g==1),Z(g==2));
    d = meanEffectSize(Z(g==1),Z(g==2),"Effect","cliff"); 
else
    p = ttest2(Z(g==1),Z(g==2));
    d = meanEffectSize(Z(g==1),Z(g==2),"Effect","cohen"); d.Effect;
end

% After correction statistics
if kstest(Z_c(g_c==1)) || kstest(Z_c(g_c==2))
    p_c = ranksum(Z_c(g_c==1),Z_c(g_c==2));
    d_c = meanEffectSize(Z_c(g_c==1),Z_c(g_c==2),"Effect","cliff");
else
    p_c = ttest2(Z_c(g==1),Z_c(g_c==2));
    d_c = meanEffectSize(Z_c(g_c==1),Z_c(g_c==2),"Effect","cohen"); d.Effect;
end

% Performance measure of clusters after correction
confusionMatrix = confusionmat(label(cutoff),idx_test(cutoff));
TP = confusionMatrix(1,1);
TN = confusionMatrix(2,2);
FP = confusionMatrix(2,1);
FN = confusionMatrix(1,2);
sensitivity_c = TP / (TP + FN);
specificity_c = TN / (TN + FP);

% Saving results to display later
res = [{'Before correction'} [sensitivity specificity] p d.Effect AUC;
 {'After correction'} [sensitivity_c specificity_c] p_c d_c.Effect AUC_c];
res = cat(1,[{''},{'Sensitivity,Specificity'},{'p-value'},{'effect'},{'AUC'}],res);
%%