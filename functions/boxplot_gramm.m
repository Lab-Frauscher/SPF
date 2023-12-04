function [p,d,gr] = boxplot_gramm(g,Z,colors,varargin)
    
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
    
    if kstest(Z(g==1)) || kstest(Z(g==2))
        p = ranksum(Z(g==1),Z(g==2));
        d = meanEffectSize(Z(g==1),Z(g==2),"Effect","cliff"); d = d.Effect;
    else
        p = ttest(Z(g==1),Z(g==2));
        d = meanEffectSize(Z(g==1),Z(g==2),"Effect","cohen"); d = d.Effect;
    end
end