function [pval,effect,gr] = boxplot_gramm(g,Z,colors,varargin)
% boxplot_gramm plots a boxplot using the gramm toolbox with the parameters 
% used to reproduce figure 8 in the manuscript. Statistical tests are
% performed as well in accordance to whether the data is normally
% distributed or not. 
%
%   [p,d,gr] = boxplot_gramm(g,Z,colors) takes as input three vectors, 
%   first of which contains the dependent variable, the second is the 
%   independent varialble, and the third are the colors for each group.
%   Additional optional inputs (2 or 3) to specify the location of the plot
%   within a subplot grid (first two inputs), and the last optional input
%   is the handle of the gramm plot
%   
%   INPUTS:     g           Nx1 vector of labels (independent variable)
%               Z           Nx1 vector of data (dependent variable)
%               colors      Mx3 vector of RGB colors corresponding the M 
%                           unique labels
%               varargin{1} Subplot row index
%               varargin{2} Subplot column index
%               varargin{3} Existing gramm plot handle
%
%   OUTPUTS:    pval        1x1 double denoting the p-value (ttest if
%                           normally distributed, ranksum if any of the 
%                           groups are not-normally distributed)
%               effect      1x1 double denoting the effect size (cohen's d 
%                           if normally distributed, cliff's d if any of the 
%                           groups are not-normally distributed)
%               gr          struct denoting the gramm plot handle
%
%   See also scatter_gramm and plotFeatureSpace.

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

    gr(i,j) = gramm('x',ones(size(Z)),'y',Z,'color',g);
    gr(i,j).set_names('y', 'Volume of SOZ resected (mm^3)');
    gr(i,j).set_text_options('interpreter', 'tex','base_size',25);
    gr(i,j).stat_boxplot('width',0.3,'notch',true);
    gr(i,j).geom_jitter('width',0.15,'dodge',0.7,'alpha',0.5);
    gr(i,j).set_color_options('map',colors,'n_color',2,'n_lightness',1);
    gr(i,j).set_point_options('base_size',10);
    gr(i,j).no_legend();
    
    % Unpaired statistical tests.
    % Computing the p-value and effect size
    if kstest(Z(g==1)) || kstest(Z(g==2))
        pval = ranksum(Z(g==1),Z(g==2));
        effect = meanEffectSize(Z(g==1),Z(g==2),"Effect","cliff"); effect = effect.Effect;
    else
        [~,pval] = ttest2(Z(g==1),Z(g==2));
        effect = meanEffectSize(Z(g==1),Z(g==2),"Effect","cohen"); effect = effect.Effect;
    end
end