% This script analyzes the result of the vSP framework which was applied
% on all patients from the MNI and CHUGA. The results are found in 
% virtual_removal_perturbation.mat  
addpath('functions');
addpath('demo_data')
addpath('perturbation')
addpath(genpath('results_figures\sp_results'))

load('virtual_removal_perturbation.mat')

% Specify center (MNI=1,CHUGA=2)
center_names = ["MNI"; "CHUGA"]; 

% center = 1; % Uncomment if running script from here
% nBoot = 1000; % Uncomment if running script from here

% Unpacking results
centers = cell2mat(data(:,end));
label = cell2mat(data(centers==center,end-1));
N1 = sum(label==1);
N2 = sum(label==2);

rho = cell2mat(data(centers==center,[3 4 5])); 
rho = [{rho(label==1,:)}; {rho(label==2,:)}];

X = [reshape(rho{1},[],1); NaN*ones(N1,1); reshape(rho{2},[],1)];
g = [ones(N1,1); 2*ones(N1,1); 3*ones(N1,1); 4*ones(N1,1); 5*ones(N2,1); 6*ones(N2,1); 7*ones(N2,1)];
MRI_cov = data(centers == center,2);
%% Figure 3: Change in spatial system after virtually removing SOZ
figure;
boxplot(X,g,'notch','off','FactorSeparator', [1])
title(center_names(center));
movegui('west')
ylim([-1 0.5]);

% Customizing boxplot line thickness
h = findobj(gca,'tag','Median');
set(h,'linestyle','-', 'linewidth',2);
set(h,'Color',[0 0 0])
h = findobj(gca,'tag','Upper Whisker');
set(h,'linestyle','-', 'linewidth',2);
h = findobj(gca,'tag','Lower Whisker');
set(h,'linestyle','-', 'linewidth',2);
h = findobj(gca,'tag','Box');
set(h,'linestyle','-', 'linewidth',2);
set(h,'Color',[0 0 0])
h = findobj(gca,'tag', 'Upper Adjacent Value');
set(h,'linestyle','-', 'linewidth',2);
set(h,'Color',[0 0 0])
h = findobj(gca,'tag', 'Lower Adjacent Value');
set(h,'linestyle','-', 'linewidth',2);
set(h,'Color',[0 0 0])

% Specifying color of boxplots to alternate
colors = [1 0 0; 0 0 1; 1 0 0; 0 0 0];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(mod(j-1,4)+1,:),'FaceAlpha',.2);
end

xticks = 1:8;
xticklabels = {'BR', 'AR', 'RR', ''};
set(gca,'xtick',xticks,'xticklabel',xticklabels,'FontSize',30)
ylabel('\rho','FontSize',30)

hold on;

N = [N1;N2];
for i = 1:2 % Looping over surgical outcomes
    z = rand(N(i),1)/10; z = z - mean(z); % Adding jitter

    hs = scatter(1 + 4*(i-1) + z, rho{i}(:,1), 36, 'MarkerFaceColor', 'r', 'jitter','off','JitterAmount',0.1, 'MarkerEdgeColor', 'k');
    hs.MarkerFaceAlpha = 0.5;

    hs = scatter(2 + 4*(i-1) + z, rho{i}(:,2), 36, 'MarkerFaceColor', 'b', 'jitter','off','JitterAmount',0.1, 'MarkerEdgeColor', 'k');
    hs.MarkerFaceAlpha = 0.5;

    hs = scatter(3 + 4*(i-1) + z, rho{i}(:,3), 36, 'MarkerFaceColor', 'r', 'jitter','off','JitterAmount',0.1, 'MarkerEdgeColor', 'k');
    hs.MarkerFaceAlpha = 0.5;    

    % Drawing lines to connect 1 and 2
    for j = 1:length(z)
        plot([ones(N(i),1) + 4*(i-1) 2*ones(N(i),1)+ 4*(i-1)] + z(j), [rho{i}(j,1) rho{i}(j,2)], 'color', [0 0 0 0.5], 'LineStyle', '--', 'LineWidth', 1);
    end    
    % Drawing lines to connect 2 and 3
    for j = 1:length(z)
        plot([2*ones(N(i),1) + 4*(i-1) 3*ones(N(i),1)+ 4*(i-1)] + z(j), [rho{i}(j,2) rho{i}(j,3)], 'color', [0 0 0 0.5], 'LineStyle', '--', 'LineWidth', 1);
    end
end

%% Statstical anlaysis 
P = []; D = [];
P = [{'SF (p-value)'} {'Non-SF (p-value)'}];
D = [{'SF (effect)'} {'Non-SF (effect)'}];
for i = 1:2 % Looping over surgical outcomes
    if kstest(rho{i}(:,1)) || kstest(rho{i}(:,2))
        [P{2,i,1},H(1,i,1)] = signrank(rho{i}(:,1), rho{i}(:,2));
        Effect = meanEffectSize(rho{i}(:,1), rho{i}(:,2), 'Effect', 'cliff'); D{2,i,1} = Effect.Effect; 
        tests{1,i,1} = 'signrank/cliff';
    else
        [~, P{2,i,1}] = ttest(rho{i}(:,1), rho{i}(:,2));
        Effect = meanEffectSize(rho{i}(:,1), rho{i}(:,2), 'Effect', 'cohen'); D{2,i,1} = Effect.Effect; 
        tests{1,i,1} = 't-test/cohen';
    end
    if kstest(rho{i}(:,1)) || kstest(rho{i}(:,3))
        [P{3,i,1},~] = signrank(rho{i}(:,1), rho{i}(:,3));
        Effect = meanEffectSize(rho{i}(:,1), rho{i}(:,3), 'Effect',  'cliff'); D{3,i,1} = Effect.Effect;  
        tests{2,i,1} = 'signrank/cliff';
    else
        [~, P{3,i,1}] = ttest(rho{i}(:,1), rho{i}(:,3));
        Effect = meanEffectSize(rho{i}(:,1), rho{i}(:,3), 'Effect',  'cohen'); D{3,i,1} = Effect.Effect; 
        tests{2,i,1} = 't-test/cohen';
    end
    if kstest(rho{i}(:,2)) || kstest(rho{i}(:,3))
        [P{4,i,1},~] = signrank(rho{i}(:,2), rho{i}(:,3));
        Effect = meanEffectSize(rho{i}(:,2), rho{i}(:,3), 'Effect',  'cliff'); D{4,i,1} = Effect.Effect;
        tests{3,i,1} = 'signrank/cliff';
    else
        [~, P{4,i,1}] = ttest(rho{i}(:,2), rho{i}(:,3));
        Effect = meanEffectSize(rho{i}(:,2), rho{i}(:,3), 'Effect',  'cohen'); D{4,i,1} = Effect.Effect; 
        tests{3,i,1} = 't-test/cohen';
    end
end

% First column is seizure-free, second is non-seizure-free
% First row -> paired test between before removal and after removal
% Second row -> paired test between before removal and random removal
% Third row -> paired test between after removal and random removal
res1 = cat(2,[{''};{'rho_BR/rho_AR'};{'rho_BR/rho_RR'};{'rho_AR/rho_RR'}],[P D]);
%% Figure 4: Surgical outcome classification using the perturbation strength
rho_br = [rho{1}(:,1); rho{2}(:,1)];
rho_ar = [rho{1}(:,2); rho{2}(:,2)];

% Computing perturbation strength
nu = 2;
rho_hat = log(nu + abs(rho_br)./(abs(rho_ar)+1e-3)); % Added 1e-3 in case denominator is zero
g = [ones(N1,1); 2*ones(N2,1)];

red =  [255 94 105]/255;
green = [157, 216, 102]/255;

% Figure 4a: Comparing perturbation strength with surgical outcome
figure;
clear gr;
gr = gramm('x', g, 'y', rho_hat, 'color', g);
gr.stat_boxplot('width',0.7,'notch',false);
gr.set_color_options('map', [green; red],'n_color',2,'n_lightness',1);
gr.set_text_options("base_size",20);
gr.set_title(center_names(center));
gr.set_names('x','','y','Perturbation strength');
gr.no_legend();
gr.axe_property('XTick',[1 2],'XTickLabel',[{'Seizure-free'}; {'Non-seizure-free'}]);
gr.draw();

% Plotting MRI positive patients with square marker
idx = strcmpi(MRI_cov,"MRI-positive");
gr.update('x', g(idx,1),'y', rho_hat(idx), 'color', g(idx,1));
gr.geom_jitter('alpha',0.5,'dodge',0.75,'width',0.3);
gr.set_point_options('base_size',10,'markers',{'s'});
gr.set_color_options('map', [green; red],'n_color',2,'n_lightness',1);
gr.no_legend();     
gr.draw();

% Plotting MRI negative patients with circle marker
idx = strcmpi(MRI_cov,"MRI-negative");
gr.update('x',g(idx,1), 'y', rho_hat(idx), 'color', g(idx,1));
gr.geom_jitter('alpha',1,'dodge',0.75,'width',0.3);
gr.set_point_options('base_size',15,'markers',{'o'});
gr.set_text_options("base_size",20);
gr.set_color_options('map', [green; red],'n_color',2,'n_lightness',1);
gr.no_legend();
gr.draw();

% Plotting threshold obtained from ROC curve
threshold_value = 1.1982; % Threshold obtained by selecting the balanced operating point on the ROC curve
gr.update('x', [0 20], 'y', [threshold_value threshold_value]);
gr.geom_line("alpha",0.5);
gr.set_color_options("map",[0 0 0],'n_color',2,'n_lightness',1);
gr.set_line_options("styles","--");
gr.axe_property('XLim',[0.5 2.5]);
gr.axe_property('YLim',[0.5 7]);
gr.draw();

set(gcf,'units','normalized','outerposition',[0 0 0.5 0.7]);
movegui('center')

% Statistical analysis
if kstest(rho_hat(g==1)) || kstest(rho_hat(g==2))
    Effect = meanEffectSize(rho_hat(g==1), rho_hat(g==2), 'Effect',  'cliff'); d = Effect.Effect;  
    p = ranksum(rho_hat(g==1),rho_hat(g==2));
    test = 'ranksum/cliff';
else
    Effect = meanEffectSize(rho_hat(g==1), rho_hat(g==2), 'Effect',  'cohen'); d = Effect.Effect;  
    p = ttest(rho_hat(g==1),rho_hat(g==2));
    test = 'ttest/cohen';
end

% Figure 4b: ROC analysis. Set center=2 at the beginning of the script to obtain Figure 4c
% Computing boostrapped ROC curve to obtain confidence intervals of AUC.
% NOTE: this can take a while!
[ROC_X,ROC_Y,ROC_T,AUC,optThreshold] = perfcurve(g,rho_hat,1,'Nboot',nBoot,'XVals',[0:0.001:1],'BootType','student');

operating_point = [0.394,0.764706];
figure;
plot(ROC_X(:,1),ROC_Y(:,1),'k-');
hold on;
% Only plot the operating point if using the training data (MNI)
if center == 1; plot(operating_point(1),operating_point(2),'rx','MarkerSize',15,'LineWidth',2); end 
% Plotting random classifier performance (diagonal)
plot(0:0.1:1,0:0.1:1,'color',[0.2 0.2 0.2 0.5],'LineStyle', '--')
ylim([0 1])

% Shading confidence intervals estimated using bootstrapping
patch([ROC_X(:,1)' fliplr(ROC_X(:,1)')], [ROC_Y(:,2)' fliplr(ROC_Y(:,3)')], [0.5, .5, .5], 'FaceAlpha',0.2, 'EdgeColor','none')
ax = gca';
ax.FontSize = 16;
ylabel('Sensitivity');
xlabel('1-Specificity');
title(center_names((center)));
movegui('east')
res2 = [{test} p d AUC];  % AUC: [mean, lower CI, upper CI]
res2 = cat(1,[{'test'},{'p-value'},{'effect'},{'AUC (mean, lower CI, upper CI'}],res2);