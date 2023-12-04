function [rho,logData,pval] = computeCorrelation(data)
% computeCorrelation  Characterizes a given spatial system as defined by
% a power law function y=ax^k as described in the original manuscript. 
%
%   [rho,logData,pval] = computeCorrelation(data) takes as input an NxL 
%   double array with N channels and L samples. The function returns rho 
%   which is the Pearson's correlation in the log-log space, and the 
%   corresponding pvalue. For visualization purposes, logData is also 
%   returned, which is the log-log transformation applied to data.
%
%   See also computeSPMap and virtualRemovalSP.
        x = data(:,1);
        y = data(:,2); 
        
        % Removing the spatial reference from correlation plot
        y(x == 0) = []; x(x == 0) = [];

        % Adding constant, as log(0) is ill-defined
        y = log(0.1+y);
        x = log(0.1+x);

        logData = [x y];

        if ~isempty(x) || ~isempty(y)
            % if computeModel; m = fitlm(x,y); end
            [rho,pval] =  corr(x,y,'Type','Pearson');
        else
            rho = 0;
            pval = Inf;
        end
        
        if isinf(rho) || isnan(rho)
            rho = 0;
        end
        
end