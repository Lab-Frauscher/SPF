% Example end-to-end pipeline of extracting data, preprocessing,
% feature extraction and spatial perturbation. The first part of the script
% reads in five minutes of MNI SEEG data that was downsampled to 200Hz due 
% to file size limitations. The data is notch filtered, and the Janca 
% detector [1] is used to detect IEDs in the SEEG signal. Then IED will be
% classified as having significant gamma activity preceding the IED using
% code previously developed in our lab [2,3]. Extracerebral and artifact 
% channels are subsequently removed before applying the virtual-removal SP
% framework. Then finally, the rSP framework is applied to the data,
% constructing the SP map. The SP map for the patient is plotted for 
% illustration purposes. 
%
% References:
% [1] Janca, R., Jezdik, P., Cmejla, R. et al. Detection of Interictal 
% Epileptiform Discharges Using Signal Envelope Distribution Modelling: 
% Application to Epileptic and Non-Epileptic Intracranial Recordings. 
% Brain Topogr 28, 172–183 (2015). https://doi.org/10.1007/s10548-014-0379-1
%
% [2] Thomas, J., Kahane, P., Abdallah, C., Avigdor, T., Zweiphenning, W.J.E.M., 
% Chabardes, S., Jaber, K., Latreille, V., Minotti, L., Hall, J., Dubeau, F., 
% Gotman, J. and Frauscher, B. (2023), A Subpopulation of Spikes Predicts 
% Successful Epilepsy Surgery Outcome. Ann Neurol, 93: 522-535. 
% https://doi.org/10.1002/ana.26548
%
% [3] https://github.com/Lab-Frauscher/Spike-Gamma
addpath('functions')
addpath("demo_data")
addpath("perturbation")
addpath("spike_gamma_code");

patient_number =1; % 1 is a seizure-free case, 2 is a non-seizure-free case
%% Step 1: Load the data of interest
% For this demo script, we will load the five minute sample data from the 
% MNI downsampled to 200Hz
demo_data = ["P4_sample.mat", "P58_sample.mat"];
load(fullfile('demo_data', demo_data(patient_number)));
load('demo_data/patient_data.mat');
patient = patient(patient_number);
fs = patient.sampling_freq;
%% Step 2: Preprocess the data
% Before running the Janca detector [1], a notch filter will be applied
% (60Hz for North American power line interference). Then a bandpass filter
% is applied before classifying the IED as an IED-gamma.
wo = 60/(fs/2);  
bw = wo/500;

[bn,an] = iirnotch(wo,bw);
[b,a] = butter(4, [0.3 70] * 2/fs, 'bandpass');

% The default settings used for the spike detector
settings = '-bl 10 -bh 60 -h 60 -jl 3.65 -dec 200';

% Zero-phase filtering 
data = filtfilt(bn,an,data);

% Applying the spike detector
out = spike_detector_hilbert_v25(data,fs,settings); % Janca detector [1]
% Post-processing detections to remove artifacts and potential spindles
out_pp = postprocessing_v2(out,fs,size(data,2))';

% Data is bandpass filtered before classifying IED as having preceding
% gamma activity
data_bp = filtfilt(b,a,data);
IED_gamma = computeSpikeGamma(data_bp,fs,out_pp,1); % Compute IED-gamma rates [2,3]

%% Step 3: Preprocessing channels
% Now that the signals are preprocessed, the channels marked as
% extracerebral, white matter or artifactual will be removed from
% subsequent analysis

% Specifying the feature vector and the coordinates for each feature
feature = IED_gamma;
coordinates = patient.MNI;

% Extracting 'bad' channel markings (i.e., extracerebral, white matter and
% artifacts)
soz_ez = patient.soz_ez;
soz = str2double(soz_ez(:,2));
invalid_indices = find(soz > 1);

% Removing invalid channels from coordinates and feature vector
feature(invalid_indices,:) = [];
coordinates(invalid_indices,:) = [];
%% Step 4: Apply ranked spatial perturbation framework 
% We can finally apply the rSP framework on the defined features and
% channel coordinates. The result of the the rSP framework is the spatial
% perturbation map, which we hypothesize can assess quality of the
% implantation of the SOZ. The function also extracts the features
% mentioned in the manuscript froms quadrants 1,2, and 4, and the
% probability model is applied to the features, giving the likelihood that
% the SOZ was well sampled.

[spMap,XBinned,yBinned,features,probability_sampled] = computeSPMap(feature,coordinates);
[{patient.name} {'Probability SOZ is well-sampled'} probability_sampled]
%% Step 5: Visualizing the results
figure;
colormap('jet');
h = pcolor(XBinned, yBinned, spMap);
set(h, 'EdgeColor', 'none');
set(gca,'FontSize',20);
title(patient.name);
xlabel("distance to perturbation centroid");
ylabel('Channel index');
colorbar;
clim([0 .5]);
