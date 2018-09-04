%% A script to generate all the figures necessary for the paper
% This file generates figures, based on pre-computed performance measures
% for each evailable dataset, illustrating statistical performance of each
% algorithm used in comparison

%% Generate plots for LASA dataset
% load aggregate results
LASA_obj = load('C:\Harish\Code\SkillLearning-using-Optimization\results\LASA_dataset\aggregate_results\LASA_aggregate_evaluation.mat');

% plot results
figure('Name','LASA Datset','units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1);
boxplot(LASA_obj.SEA_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(LASA_obj.SEA_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Swept Error Area (m^2)')
title('Swept Error Area (m^2)','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

subplot(1,3,2);
boxplot(LASA_obj.SSE_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(LASA_obj.SSE_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Sum of Squared Distances (m^2)')
title('Sum of Squared Distances (m^2)','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

subplot(1,3,3);
boxplot(LASA_obj.DTWD_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(LASA_obj.DTWD_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Dynamic Time Warping Distance')
title('Dynamic Time Warping Distance','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

%% Generate plots for RAIL picking datset
% load aggregate results
RAIL_picking_obj = load('C:\Harish\Code\SkillLearning-using-Optimization\results\RAIL_dataset\picking\subject_5_picking_evaluation.mat');

% plot results
figure('Name','RAIL picking Datset','units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1);
boxplot(RAIL_picking_obj.SSE_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(RAIL_picking_obj.SSE_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Sum of Squared Distances (m^2)')
title('Sum of Squared Distances (m^2)','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

subplot(1,2,2);
boxplot(RAIL_picking_obj.DTWD_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(RAIL_picking_obj.DTWD_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Dynamic Time Warping Distance')
title('Dynamic Time Warping Distance','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

%% Generate plots for RAIL pressing datset
% load aggregate results
RAIL_pressing_obj = load('C:\Harish\Code\SkillLearning-using-Optimization\results\RAIL_dataset\pressing\subject_5_pressing_evaluation.mat');

% plot results
figure('Name','RAIL Pressing Datset','units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1);
boxplot(RAIL_pressing_obj.SSE_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(RAIL_pressing_obj.SSE_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Sum of Squared Distances (m^2)')
title('Sum of Squared Distances (m^2)','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

subplot(1,2,2);
boxplot(RAIL_pressing_obj.DTWD_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(RAIL_pressing_obj.DTWD_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Dynamic Time Warping Distance')
title('Dynamic Time Warping Distance','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

%% Generate plots for RAIL pushing datset
% load aggregate results
RAIL_pushing_obj = load('C:\Harish\Code\SkillLearning-using-Optimization\results\RAIL_dataset\pushing\subject_5_pushing_evaluation.mat');

% plot results
figure('Name','RAIL Pushing Datset','units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1);
boxplot(RAIL_pushing_obj.SSE_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(RAIL_pushing_obj.SSE_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Sum of Squared Distances (m^2)')
title('Sum of Squared Distances (m^2)','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)

subplot(1,2,2);
boxplot(RAIL_pushing_obj.DTWD_list,'symbol','')
set(findobj(gca,'type','line'),'linew',2)
hold on
plot(mean(RAIL_pushing_obj.DTWD_list),'p','color',[0 0.5 0],'LineWidth',3,'MarkerSize',10)
set(gca,'FontSize',25)
% ylabel('Dynamic Time Warping Distance')
title('Dynamic Time Warping Distance','FontSize',20)
xticklabels({'Position','Tangent','Laplace','Uniform','Learned'})
xtickangle(45)
