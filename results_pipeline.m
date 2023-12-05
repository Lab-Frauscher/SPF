% The purpose of this script is to reproduce all the results in the
% manuscript. Just click run and enjoy! 
% Note: it's more preferable to run section by section, or your screen will
% be flooded. Don't say I didn't warn you :)

% You can directly run the scripts individually instead of using this 
% automated pipeline if you wsh to play with some parameters, have fun!

clear all;
clc;
close all;

addpath(genpath('results_figures'));
%% Figures 3 and 4: Illustrates the change in the spatial system after 
% virtually removing the SOZ for seizure-free patients (Figure 3a) and
% non-seizure-free patients (Figure 3b). The perturbation strength was
% calculated and shown to be significantly larger in seizure-free patients
% for both centers (Figure 4)

% I selected 100 to speed things up, 1000 was used in the manuscript
nBoot = 100; %
results_1 = [];
results_2 = [];
for center = 1:2
    figure3_4_perturbation_analysis;
    results_1 = [results_1; res1];
    results_2 = [results_2; res2];
end
% Figure 3
disp('Figure 3a (MNI)')
disp(results_1(1:4,:));
disp('Figure 3b (CHUGA)')
disp(results_1(5:end,:));

% Figure 4 
disp('Figure 4a')
results_2 = cat(2,[{''};{'MNI'};{'CHUGA'}],[results_2(1:2,:); results_2(end,:)]);
disp(results_2)
%% Figure 6: Illustrates the spatial perturbation map of two patients, one 
% seizure-free (P4) and the other is non-seizure-free (P58)
patient_number=1; % Seizure-free (P4)
figure6_sp_map_examples; movegui('center');

patient_number=2; % Non-seizure-free (P58)
figure6_sp_map_examples; movegui('east')

%% Figure 7: Illustrates the clustering analysis for MNI and Chuga as well 
% as the results of the probability model on all patients, and after
% correcting for patients with incomplete resections of the SOZ
center = 1;
figure7_cluster_analysis; movegui('north'); res_1 = res;
center = 2;
figure7_cluster_analysis; movegui('south'); res_2 = res;

disp('Figure 7a,c,e (MNI)')
disp(res_1); % MNI results
disp('Figure 7b,d,f (CHUGA)')
disp(res_2); % CHUGA results

%% Figure 8 and 9: Illustrates the differences in SOZ volume when comparing
% surgical outcome within the good sampling cluster, and when comparing 
% across clusters (Figure 8). The SOZ volume was compared after correcting 
% for patients with insufficient resections of the SOZ in both centers.
% Figure 9 demonstrates the inverse correlation between percent SOZ 
% removed and the distance from the patients to the good sampling centroid
% in the feature space (in seizure-free patients only). 
figure8_9_covariate_analysis;
disp(cat(1,[{''},{'p-value'},{'Spearman''s rho'}],results))