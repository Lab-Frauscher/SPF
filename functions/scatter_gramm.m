function [pval,rho,gr] = scatter_gramm(g,Z,labels,colors,color_fp,fp,varargin)
% scatter_gramm plots a scatter plot using the gramm toolbox with the 
% parameters used to reproduce figure 9 in the manuscript. Statistical 
% tests are performed as well in accordance to whether the data is normally
% distributed or not. 
%
%   [pval,rho,gr] = scatter_gramm(g,Z,colors) takes as input four vectors, 
%   first of which contains the dependent variable, the second is the 
%   independent variable, the third is a vector of labels for each point
%   in the scatter plot, and the fourth vector is an RGB triplet for each 
%   label (Mx3 for M unique number of labels). Additional optional inputs 
%   (2 or 3) to specify the location of the plot within a subplot grid 
%   (first two inputs), and the last optional input is the handle of the 
%   gramm plot
%   
%   INPUTS:     g           Nx1 independent variable vector
%               Z           Nx1 dependent variable vector
%               labels      Nx1 vector of labels (integer or boolean)
%               colors      Mx3 vector of RGB colors corresponding the M 
%                           unique labels
%               varargin{1} Subplot row index
%               varargin{2} Subplot column index
%               varargin{3} Existing gramm plot handle
%
%   OUTPUTS:    pval        1x1 double denoting the p-value (ttest if
%                           normally distributed, ranksum if any of the 
%                           groups are not-normally distributed)
%               rho         1x1 double denoting the Spearman's rank
%                           correlation coefficient
%               gr          struct denoting the gramm plot handle
%
%   See also boxplot_gramm and plotFeatureSpace.
    
    if length(varargin) < 3
        i = varargin{1};
        j =varargin{2};
    elseif length(varargin) == 3
        i = varargin{1};
        j =varargin{2};
        gr = varargin{3};
    else
        error('Too many inputs!');
    end

    gr(i,j) = gramm('x',g,'y',Z,'color',labels);
    gr(i,j).set_color_options('map',colors,'n_color',2,'n_lightness',1);
    gr(i,j).geom_point();
    gr(i,j).set_point_options('base_size',10);
    gr(i,j).set_text_options('interpreter', 'tex','base_size',20);
    gr(i,j).no_legend();
    gr(i,j).draw();

    gr(i,j).update('x',g(labels==fp),'y',Z(labels==fp));
    gr(i,j).set_color_options('map',color_fp,'n_color',1,'n_lightness',1);
    gr(i,j).geom_point();
    gr(i,j).set_point_options('base_size',5);
    gr(i,j).no_legend();
    
    % Testing for correlations in the scatter plot
    [rho,pval]=corr(g,Z,'type','Spearman');
end