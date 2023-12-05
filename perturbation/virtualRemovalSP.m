function [rhos,pvals,spatial_system,flag] = virtualRemovalSP(feature,coordinates,soz)
% computeSPMap  Applies a perturbation on the spatial system by virtually
% removing the SOZ as described in the original manuscript. 
%
%   [rhos,pvals,data] = virtualRemovalSP(feature,coordinates,soz) takes as
%   input the channel-level feature values (Nx1 double array) for a segment
%   of interest, the coordinates of each channel (Nx3 double array), and
%   the variable soz (boolean array) which indicates whether a channel is
%   is in the SOZ.
%   virtualRemovalSP returns the correlation values in 'rhos' as a 3x1 
%   double array, which chracterizes the spatial system before removal of 
%   the SOZ, after removal of the SOZ, and after random removal of 
%   non-SOZ channels. p-values for the correlation are also returns in the 
%   second output argument as a 3x1 double array. For visualization 
%   purposes, the function also returns the spatial system before after 
%   virtually removing the SOZ as a 2x2 cell array. Each cell is an Nx2 
%   double array  
%   - spatial_system{1,1} contains distances and feature values before 
%   removal of the SOZ
%   - spatial_system{1,2} is the spatial system in the log-log space
%   (before removal)
%   - spatial_system{2,1} contains distances and feature values after 
%   removal of the SOZ
%   - spatial_system{2,2} is the spatial system in the log-log space
%   (after removal)
%
%   INPUTS:     feature             Nx1 vector  (N channels)
%               coordinates         Nx3 vector  (three dimensions)
%               soz                 Nx1 boolean vector indicating whether a given channel belongs to the SOZ        
%   OUTPUTS:    rhos                3x1 vector of Pearson's correlations applied to the spatial system
%                                       - before removal of the SOZ
%                                       - after removal of the SOZ
%                                       - random removal of non-SOZ channels
%               pvals               3x1 vector of p-values corresponding to rho
%               spatial_system      2x2 cell array of spatial systems before and after virtually removing the SOZ
%               flag                1x1 boolean indicating whether the segment should be used
%
%   See also computeCorrelation and computeSPMap.

% Initializing outputs
rhos = []; pvals = []; logData = {[]}; flag = false;
%% Step 1: Constructing and characterizing spatial system
% The spatial system is constructing by coupling channel-level features
% with their distances to a spatial reference (defined as the channel with 
% the highest feature value in the SOZ). The spatial system is
% characterized by the goodness-of-fit of the power-law function using the 
% MATLAB function 'computeCorrelation.m'

sozIdx = find(soz==1);
[~,ii] = max(feature(sozIdx,:),[],1);  % Computing spatial reference
varphi_sr_soz = sozIdx(ii)';

x = pdist2(coordinates,coordinates(varphi_sr_soz,:));
y = feature;
spatial_system{1,1} = [x y]; % For visualizing results later

% Checking max rates for the given segment
rates = max(feature,[],1);
if rates > 1 % Only continue if there are sufficient IED-gammas
    % Characterizing spatial system before virtual removal (BR) of SOZ
    [rho_br,logData,pval_br]=computeCorrelation([x,y]);
    spatial_system{1,2} = logData; % For visualizing results later
    %% Step 2: Applying perturbation by virtually removing the SOZ
    % The SOZ is virtually removed, and the spatial reference is re-computed
    % relative to the channel with the highest feature value (outside the
    % SOZ). The perturbed spatial system is characterized as done in Step 1

    % Virtual removal of the SOZ
    feature_nez = feature(soz==0,:); 
    coordinates_nsoz = coordinates(soz==0,:);
    
    % Computing spatial reference
    [~,varphi_sr_nsoz] = max(feature_nez(sozIdx,:),[],1);  
    
    x = pdist2(coordinates_nsoz,coordinates_nsoz(varphi_sr_nsoz,:));
    y = feature_nez;

    % Characterizing spatial system after virtual removal (AR) of SOZ
    [rho_ar,logData,pval_ar]=computeCorrelation([x,y]);
    spatial_system{2,1} = [x y]; % For visualizing results later
    spatial_system{2,2} = logData; % For visualizing results later

    % Checking max rates for each segment
    rates = max(feature_nez,[],1);
    % Setting correlation to zero if segment has low rates
    rho_ar = rho_ar * (rates > 1);
    
    %% Step 3: Randomly removing non-SOZ channels to compare with spatial 
    % system before virtually removing the SOZ

    N = 100;        % Values used in the manuscript
    ratio = 0.2;    % Values used in the manuscript
    rng('default'); % For reproducability
    
    idxToSelect = (1:size(feature,1))';
    % Only randomly select from non-SOZ indices 
    idxToSelect = setdiff(idxToSelect, sozIdx); 
    L = round(length(sozIdx)*ratio);
    
    rho_rr = [];
    pval_rr = [];
    for i = 1:N % Looping over all bootstraps
        % Obtain random indices of non-SOZ channels
        rand_indices = setdiff(1:size(feature,1), randsample(idxToSelect,L)');

        % Removing random non-SOZ channels
        feature_rr = feature(rand_indices,:); 
        coordinates_rr = coordinates(rand_indices,:);
    
        % Computing spatial reference
        [~,ii] = max(feature_rr(sozIdx,:),[],1);  
        varphi_sr_rr = sozIdx(ii)';
        
        x = pdist2(coordinates_rr,coordinates_rr(varphi_sr_rr,:));
        y = feature_rr;
        
        % Checking max rates for each segment
        rates = max(feature,[],1);
        lowRates = rates <= 1;
    
        % Characterizing spatial system after random removal of non-SOZ
        % channels
        [rho_rr_tmp,~,pval_rr_tmp]=computeCorrelation([x,y]);

        rho_rr = [rho_rr; rho_rr_tmp(~lowRates)]; % Only consider segments with rates > 1
        pval_rr = [pval_rr; pval_rr_tmp(~lowRates)];
    end
    %Computing median over bootsraps
    rho_rr = mean(rho_rr);
    pval_rr = mean(pval_rr);

    % Saving results to output
    rhos = [rho_br;rho_ar;rho_rr];
    pvals = [pval_br;pval_ar;pval_rr];
else
    % Flag set to 1 to indicate that the segment should not be included
    flag = true;
end