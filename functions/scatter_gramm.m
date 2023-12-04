function [pval,rho,gr] = scatter_gramm(g,Z,labels,colors,color_fp,fp,varargin)
    
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
    
    [rho,pval]=corr(g,Z,'type','Spearman');
end