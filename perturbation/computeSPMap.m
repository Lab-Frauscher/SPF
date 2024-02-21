function [spMap,XBinned,yBinned,spMap_features,probability_sampled] = computeSPMap(feature,coordinates)
% computeSPMap  Applies a series of perturbations to the spatial system
% described in the original manuscript. Perturbation strengths are
% spatially ranked to produce a spatial perturbation (SP) map. Features are
% computed and a classification is made on the implantation scheme
% (well-sampled or poorly sampled)
%
%   [spMap,XBinned,yBinned,features,probability_sampled] =
%   computeSPMap(feature,coordinates) produces the SP map using the
%   channel-level features (NxM) and their corresponding coordinates (Nx3),
%   where N is the number of channels, and M is the number of
%   observations. The function returns spMap, an Nx200 matrix of
%   spatially ranked perturbation strengths (the SP map). The x and y axes 
%   are returned as XBinned (200x1) and YBinned(Nx1). The variable
%   'features' (2x2) is extracted from the SP map
%       features(1,1) = mean positive perturbation indices of quadrant 1
%       features(1,2) = mean positive perturbation indices of quadrant 2
%       features(2,1) = mean positive perturbation indices of quadrant 3
%       features(2,2) = mean positive perturbation indices of quadrant 4
%   The features are used to assign a probabillity using the model
%   developed (see figure7_cluster_analysis.m) as described in the
%   manuscript. This is returned as 'probability_sampled' (1x1 double)
%
%   INPUTS:     feature             NxM vector  (N channels, M observations)
%               coordinates         Nx3 vector  (N channels, three dimensions)
%   OUTPUTS:    spMap               Nx200 vector of spatially ranked perturbation strengths
%               XBinned             200x1 vector of binned distances
%               YBinned             Nx1 vector of channel indices (channel rankings)
%               spMap_features      2x2 matrix of features computed on quadrants 1,2,3 and 4 
%               probability_sampled 1x1 double ranging from 0 to 1, indicating probability of sampling the SOZ 
%
%   See also computeCorrelation and virtualRemovalSP.

    % Applying the method to each segment
    pStrength = zeros(size(feature,1),size(feature,2)); 
    for seg = 1:size(feature,2) % Computing perturbation strength for each segment
        feature_values = feature(:,seg);

        %% Step 1: Construct spatial system
        % Compute spatial reference
        [~,varphi_sr] = max(feature_values); 
 
        dist = pdist2(coordinates,coordinates(varphi_sr,:));
        rho_ref = computeCorrelation([dist feature_values]);
        pStrength(varphi_sr,seg) = 1;

        % Remove initial spatial reference from further analysis 
        % (we want to apply series of perturbation and evaluate change from
        % initial spatial system
        idxs = setdiff(1:size(feature,1),varphi_sr); 

        %% Step 2: Permute spatial reference
        for i = idxs
            varphi_sr_i = i; % Permuting spatial reference            
            
            x = pdist2(coordinates(idxs,:),coordinates(varphi_sr_i,:));
            y = feature_values(idxs);

            % Characterizing perturbed spatial system
            [rho,~,~] = computeCorrelation([x y]); 
            pStrength(i,seg) = 1-(rho-rho_ref);
        end
    end
    % Computing median over all segments
    pStrength = abs(median(pStrength,2));

    % Determine region of high perturbation strength
    high = prctile(pStrength,70); 
    coord_roi = coordinates(pStrength>=high,:);
   
    % Compute centroid of region with high perturbation strength
    varphi_ps = mean(coord_roi,1);
    
    %% Step 3: Spatial ranking (second-step perturbation)
    % Ranking distances of each channel to region with high perturbation
    distances = pdist2(coordinates,varphi_ps);
    % Define the number of bins for each coordinate dimension
    numBinsX = 200;  % Number of bins along the X-axis
    % Compute the bin indices for each observation
    binIndicesX = discretize(distances, linspace(min(distances(:, 1)), max(distances(:, 1)), numBinsX + 1));

    XBinned = (min(binIndicesX):max(binIndicesX))';
    yBinned = (1:size(feature,1))';
    spMap = zeros(size(feature,1),length(XBinned));

    % Binning data, and averaging perturbation strengths coinciding with the same bin 
    for j = 1:length(XBinned)
        idx = find(binIndicesX == j);
        if ~isempty(idx)
            spMap(idx,j) = ( spMap(idx,j) + pStrength(idx) )/2;
        end
    end

    % Ranking the channels by their perturbation strength
    [~,isort] = sort(pStrength);
    spMap = spMap(isort,:);

    % Applying morphological operations to close points in proximity
    for angle = 10:0.5:80
        spMap = imclose(spMap,strel('line', 20, angle));
    end
    for angle = 80:-0.5:10
        spMap = imclose(spMap,strel('line', 20, angle));
    end
    
    %% Step 4: Extracting features from the SP map
    % Features are extracted from the SP map to quantify the implantation
    % scheme
    spMap_features = [];
    xMax = max(XBinned);
    yMax = size(spMap,1);
    
    Q1 = spMap(yBinned >= yMax/2, XBinned <= xMax/2); 
    spMap_features(1,1) = mean(Q1(Q1 > 0)); 

    Q2 = spMap(yBinned >= yMax/2, XBinned > xMax/2); 
    spMap_features(1,2) = mean(Q2(Q2 > 0)); 

    Q3 = spMap(yBinned < yMax/2, XBinned <= xMax/2); 
    spMap_features(2,1) = mean(Q3(Q3 > 0)); 

    Q4 = spMap(yBinned < yMax/2, XBinned > xMax/2); 
    spMap_features(2,2) = mean(Q4(Q4 > 0)); 

    spMap_features(isnan(spMap_features)) = 0;
    %% Step 5: Applying probability model to SP map feature space
    % Lastly, the probability model is applied to the extracted features 
    % indicating how well the SOZ was implanted by the SEEG

    % Centroids that were estimated from the training data (MNI)
    centroids = [0.3319    0.1137    0.0696;
    0.4401    0.3843    0.2468];

    % Parameters of the model estimated from training data (MNI)
    mu_1 =  0.1202;
    mu_2 = 0.1242;
    sigma_1 = 0.0485;
    sigma_2 =  0.0604;
    v_1 = 23;
    v_2 = 25;

    % Computing t-scores
    tdist2T = @(t,v) (1-betainc(v./(v+t.^2),v/2,0.5)); % CDF of Student's T distribution

    X = [spMap_features(1,1) spMap_features(1,2) spMap_features(2,2)];
    distances_c1 = pdist2(centroids(1,:),X); % Distance to good sampling centroid
    distances_c2 = pdist2(centroids(2,:),X); % Distance to poor sampling centroid

    t_score_c1 = (distances_c1 - mu_1)/(sigma_1); % Normalized t score
    t_score_c2 = (distances_c2 - mu_2)/(sigma_2); % Normalized t score

    probability_sampled_c1 = tdist2T(t_score_c1,v_1); % Probability patient is not in good sampling centroid (c1 p-value)
    probability_sampled_c2 = tdist2T(t_score_c2,v_2); % Probability patient is not in poor sampling centroid (c2 p-value)

    probability_sampled = 1/2.*(1-probability_sampled_c1) + 1/2.*probability_sampled_c2; % Averaged probability
end