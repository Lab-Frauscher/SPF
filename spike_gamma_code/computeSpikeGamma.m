function [gammaSpikeRate,gammaSpikeFeature] = computeSpikeGamma(data,fs,spikes,verbose)
    if verbose; f = waitbar(0, "Processing data"); pause(0.5); end

    nChannels = size(data,2);
    for chIndex = 1:nChannels
        if verbose
            waitbar(chIndex/nChannels, f, "computing spike-gamma [" + string(chIndex) + "/" + nChannels + "]");
        end
        ch = data(:,chIndex);
        spike_ch = round(spikes{chIndex}*fs);
    
        % [-75 ms 225 ms]
        onset = 75e-3 * fs;
        offset = 225e-3 * fs;
        
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
        
        %% Compute gamma 
        for spIdx = 1:length(spike_ch)
            onset = 1 * fs;
            offset = 1 * fs;
            
            spike = ch(spikePeaks(spIdx) - onset : spikePeaks(spIdx) + offset);
            
            % spike_segment = ch(spikePeaks(1) - onset: spikePeaks(1) + offset);
            % spikeGamma = filtfilt(b_gamma,a_gamma,spike_segment);
            
            ref = spikePeaks - onset;
            P1 = spikeOnsets - ref;
            N2 = spikeOffsets - ref;
            
            gammaSpikeFeature{chIndex,1}(spIdx,:) = compute_gamma(spike,fs, P1(spIdx), N2(spIdx));
        end
        if ~ismissing(spike_ch)
            tmp = gammaSpikeFeature{chIndex,1} > 0;
            gammaSpikeRate(chIndex,1) = sum(tmp(:,1)) / (size(data,1)/(fs*60));
        else
            gammaSpikeRate(chIndex,1) = 0;
        end
    end
    close(f)
end