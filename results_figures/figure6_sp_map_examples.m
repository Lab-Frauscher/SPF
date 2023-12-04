addpath('functions');
addpath('demo_data')
addpath('perturbation')
addpath(genpath('results_figures\sp_results'))

load('patient_data.mat');

% 1 is a seizure-free case, 2 is a non-seizure-free case
% patient_number = 1 % Uncomment this line if running script from here
patient = patient(patient_number);
%% Unpacking spike gamma rates computed from MNI dataset
feature = patient.spike_gamma_rates; 
coordinates = patient.MNI;
%% Removing bad channels
soz_ez = patient.soz_ez;
soz = str2double(soz_ez(:,2));
invalid_indices = find(soz > 1);

% Removing invalid channels
feature(invalid_indices,:) = [];
soz_ez(invalid_indices,:) = []; 
soz(invalid_indices) = [];
coordinates(invalid_indices,:) = [];

%% Computing SP map
[spMap,XBinned,yBinned,features,probability_sampled] = computeSPMap(feature,coordinates);
disp([{patient.name} {'Probability SOZ is well-sampled'} probability_sampled]);

%% Figure 6: Plotting SP map
figure;
colormap('jet');
h = pcolor(XBinned, yBinned, spMap);
set(h, 'EdgeColor', 'none');
colorbar;
clim([0 .5]);
title(patient.name)
set(gca,'FontSize',20)