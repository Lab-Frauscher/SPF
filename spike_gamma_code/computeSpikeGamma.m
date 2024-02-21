function [gammaSpikeRate,gammaSpikeStats] = computeSpikeGamma(data,fs,spikes,verbose)
% computeSpikeGamma computed the spike-gamma rates over all channels of a 
% multichannel data using the spikes detected from the Janca detector [1].
% The spike-gamma code is implemented as described in [2,3].
%
%   [gammaSpikeRate,gammaSpikeFeature] = computeSpikeGamma(X,label,centroid) 
%   takes as input an LxN multichannel data, the sampling rate, and the 
%   spike detections from the Janca detector [1]
%   
%   INPUTS:     data            LxN data vector (L time samples, N channels)
%               fs              1x1 double denoting sampling rate (samples/sec)
%               spikes          1xN cell array containing spike-peak locations
%                               (output of postprocessing_v2.m)
%               verbose         1x1 boolean on verbosity of code execution
%
%   OUTPUTS:    gammaSpikeRate  Nx1 vector of spike-gamma rates for each
%                               channel 
%               gammaSpikeStats Nx1 cell array of statistics, each cell 
%                               contains A 3-element array containing the 
%                               maximum gamma power, the gamma frequency
%                               corresponding to the maximum power, and the 
%                               duration of the gamma activity in
%                               milliseconds [3]. If no gamma activity is 
%                               detected, it returns [0, 0, 0].
%
%   See also compute_gamma, compute_spike_boundary, spike_detector_hilbert_v25, postprocessing_v2.
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

    % Setting up progress bar
    if verbose; f = waitbar(0, "Processing data"); pause(0.5); end

    nChannels = size(data,2);
    % Looping over each channel
    for chIndex = 1:nChannels
        if verbose
            waitbar(chIndex/nChannels, f, "computing spike-gamma [" + string(chIndex) + "/" + nChannels + "]");
        end
        ch = data(:,chIndex);
        spike_ch = round(spikes{chIndex}*fs);
    
        % [-75 ms 225 ms]
        onset = round(75e-3 * fs);
        offset = round(225e-3 * fs);
        
        spike_onset = spike_ch - onset;
        spike_offset = spike_ch + offset;
        
        spikePeaks = [];
        spikeOnsets = [];
        spikeOffsets = [];
        toRemove = [];
        for j = 1:length(spike_ch)
            spike_segment = ch(spike_onset(j) : spike_offset(j));  
            [p1, n1, n2]=compute_spike_boundary(spike_segment,fs);
    
            if ~isempty(p1) % Ignore spikes where the onset wasn't detected
                spikePeaks = [spikePeaks; n1];
                spikeOnsets = [spikeOnsets; p1];
                if ~isempty(n2) % Fix the spike offset to 150 ms away from peak
                    spikeOffsets = [spikeOffsets; n2];
                else
                    spikeOffsets = [spikeOffsets; n1+round(150e-3*fs)];
                end
            else
                toRemove = [toRemove; j];
            end
        end
        spike_ch(toRemove) = [];
        spikePeaks = spikePeaks + spike_ch' - onset;
        spikeOnsets = spikeOnsets + spike_ch' - onset;
        spikeOffsets = spikeOffsets + spike_ch' - onset;
        
        %% Checking for significant gamma activity preceding spike onset 
        for spIdx = 1:length(spike_ch)
            onset = 1 * fs;
            offset = 1 * fs;
            
            spike = ch(spikePeaks(spIdx) - onset : spikePeaks(spIdx) + offset);
            
            ref = spikePeaks - onset;
            P1 = round(spikeOnsets - ref);
            N2 = round(spikeOffsets - ref);
            
            gammaSpikeStats{chIndex,1}(spIdx,:) = compute_gamma(spike,fs, P1(spIdx), N2(spIdx));
        end
        if ~ismissing(spike_ch)
            tmp = gammaSpikeStats{chIndex,1} > 0;
            gammaSpikeRate(chIndex,1) = sum(tmp(:,1)) / (size(data,1)/(fs*60));
        else
            gammaSpikeRate(chIndex,1) = 0;
        end
    end
    close(f)
end