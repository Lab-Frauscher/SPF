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
%% Figures 3 and 4

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
%% Figure 6
patient_number=1; % Seizure-free (P4)
figure6_sp_map_examples; movegui('center');

patient_number=2; % Non-seizure-free (P58)
figure6_sp_map_examples; movegui('east')

%% Figure 7
center = 1;
figure7_cluster_analysis; movegui('north'); res_1 = res;
center = 2;
figure7_cluster_analysis; movegui('south'); res_2 = res;

disp('Figure 7a,c,e (MNI)')
disp(res_1); % MNI results
disp('Figure 7b,d,f (CHUGA)')
disp(res_2); % CHUGA results

%% Figure 8 and 9
figure8_9_covariate_analysis;
disp(cat(1,[{''},{'p-value'},{'Spearman''s rho'}],results))