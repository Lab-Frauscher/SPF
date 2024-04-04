% Example end-to-end pipeline of extracting data, preprocessing,
% feature extraction and spatial perturbation. The first part of the script
% reads in five minutes of MNI SEEG data that was downsampled to 200Hz due 
% to file size limitations. The data is notch filtered, and the Janca 
% detector [1] is used to detect IEDs in the SEEG signal. Then IED will be
% classified as having significant gamma activity preceding the IED using
% code previously developed in our lab [2,3]. Extracerebral and artifact 
% channels are subsequently removed before applying the virtual-removal SP
% framework. Then finally, the scatter plots are visualized to demonstrate
% the perturbation of the spatial system after virtually removing the SOZ.
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

patient_number = 1; % 1 is a seizure-free case, 2 is a non-seizure-free case
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
feature = IED_gamma;    % Can replace with patient.spike_gamma_rates for results from all patients
coordinates = patient.MNI;

% Extracting 'bad' channel markings (i.e., extracerebral, white matter and
% artifacts)
soz_ez = patient.soz_ez;
soz = str2double(soz_ez(:,2));
invalid_indices = find(soz > 1);

% Removing invalid channels from coordinates and feature vector
feature(invalid_indices,:) = [];
soz_ez(invalid_indices,:) = []; 
soz(invalid_indices) = [];
coordinates(invalid_indices,:) = [];

%% Step 4: The virtual-removal spatial perturbation framework
% We can finally apply the vSP framework on the defined features and
% channel coordinates

[rhos,pvals,data] = virtualRemovalSP(feature,coordinates,soz);

rho_br = rhos(1);
rho_ar = rhos(2);
rho_rr = rhos(3);

% Displaying the spatial system before and after virtually removing the SOZ
disp([[{'rho_BR'}, {'rho_AR'}, {'rho_RR'}]; num2cell([rho_br rho_ar rho_rr])])

% Computing perturbation strength
pStrength = log(2 + abs(rho_br)/(abs(rho_ar)+1e-3)); % The stronger the perturbation, the better the implantation
disp([{'Perturbation strength'} pStrength]);

%% Step 5: Visualizing results
figure;
clear gr;
gr(1,1) = gramm('x',data{1,1}(:,1), 'y', data{1,1}(:,2), 'color', soz==0);
gr(1,1).geom_point('alpha',0.7);
gr(1,1).set_point_options("base_size",10);
gr(1,1).set_names('x', 'distance to spatial ref (mm)', 'y', "IED-\gamma rates min^{-1}");
gr(1,1).set_text_options('interpreter', 'tex','base_size',30);
gr(1,1).no_legend();

[~,ii]=max(feature(soz==1));
sozIdx = find(soz==1);
soz(sozIdx(ii))=[];
gr(2,1) = gramm('x', data{1,2}(:,1), 'y', data{1,2}(:,2), 'color', soz==0);
gr(2,1).geom_point('alpha',0.7);
gr(2,1).set_point_options("base_size",10);
gr(2,1).set_names('x', 'ln(distance to spatial ref)', 'y', "ln(IED-\gamma rates min^{-1})");
gr(2,1).set_text_options('interpreter', 'tex','base_size',30);
gr(2,1).no_legend();

gr(1,2) = gramm('x', data{2,1}(:,1), 'y', data{2,1}(:,2));
gr(1,2).geom_point('alpha',0.7);
gr(1,2).set_point_options("base_size",10);
gr(1,2).set_names('x', 'distance to spatial ref (mm)', 'y', "IED-\gamma rates min^{-1}");
gr(1,2).set_text_options('interpreter', 'tex','base_size',30);
gr(1,2).set_color_options('map',[76 208 225]/255,'n_color',1,'n_lightness',1);
gr(1,2).no_legend();

gr(2,2) = gramm('x', data{2,2}(:,1), 'y', data{2,2}(:,2));
gr(2,2).geom_point('alpha',0.7);
gr(2,2).set_point_options("base_size",10);
gr(2,2).set_names('x', 'ln(distance to spatial ref)', 'y', "ln(IED-\gamma rates)");
gr(2,2).set_text_options('interpreter', 'tex','base_size',30);
gr(2,2).set_color_options('map',[76 208 225]/255,'n_color',1,'n_lightness',1);
gr(2,2).no_legend();

gr.draw();
set(gcf,'Position',[559,73,1115,877]);