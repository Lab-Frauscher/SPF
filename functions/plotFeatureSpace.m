function g = plotFeatureSpace(X,label,centroid)
% plotFeatureSpace produces a scatter plot of the feature space with the 
% parameters used to figure 7 in the manuscript. 
%
%   g = plotFeatureSpace(X,label,centroid) takes as input
%   an Nx3 vector containing the features and an Nx1 vector containing the
%   labels (ie., surgical outcome) for each observation (i.e., patient). 
%   The centroid estimated using clustering techniques is also taken as input 
%   for visualization purposes. The function returns the gramm plot handle
%   
%   INPUTS:     X           Nx1 feature vector (N observations)
%               label       Nx1 feature vector (N observations)
%               centroid    Kx3 centroid vector (K clusters)
%
%   OUTPUTS:    g           struct containing the gramm plot handle
%
%   See also boxplot_gramm and scatter_gram.

    figure;
    clear g
    
    [~,idx_test] = pdist2(centroid,X,'euclidean','smallest',1); idx_test = idx_test'; 

    factor = 2;
    
    labelsString = string(idx_test);
    labelsString (idx_test == 1) = "Cluster 1";
    labelsString (idx_test == 2) = "Cluster 2";
    
    red =  [255 94 105]/255;
    green = [157, 216, 102]/255;
    
    g=gramm('x',X(:,1),'y',X(:,2),'z',X(:,3),'color',cellstr(labelsString(:)));
    g.set_names('x','Q1','y','Q2','z','Q4');
    g.geom_point("alpha",1);
    g.set_text_options('interpreter', 'tex','base_size',20*factor);
    g.set_color_options('map',[green;red],'n_color',2,'n_lightness',1);
    g.set_point_options("base_size",10*factor);
    g.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5]);
    g.no_legend();

    g.update('x',centroid(:,1),'y',centroid(:,2),'z',centroid(:,3));
    g.set_names('x','Q1','y','Q2','z','Q4');
    g.geom_point("alpha",1);
    g.set_point_options('markers',{'s'},"base_size",20*factor);
    g.set_color_options('map',[0 0 0],'n_color',1,'n_lightness',1);
    g.no_legend();


    idxFP = label==2 & idx_test==1;
    idxFN = label==1 & idx_test==2;

    Xtmp = [X(idxFP,:); X(idxFN,:)];
    labels = [ones(sum(idxFP),1); 2*ones(sum(idxFN),1)];
    g.update('x',Xtmp(:,1),'y',Xtmp(:,2),'z',Xtmp(:,3),'color',labels);
    g.set_names('x','Q1','y','Q2','z','Q4');
    g.geom_point("alpha",1);
    g.set_point_options('markers',{'o'},"base_size",5*factor);
    g.set_color_options('map',[red;green],'n_color',2,'n_lightness',1);
    g.no_legend();
  
    axis equal
end