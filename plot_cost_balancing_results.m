%% plot cost balancing results comparing to vanilla GMM
% this function compares the results from the cost balancing for lasa
clc, clear, close all
%% load the data from file
skillnum = 15;
filenamesaved = ['skill_' num2str(skillnum) '_trained.mat'];
load(filenamesaved); % will include the followings : 'Demos','Gmms','Sols','w1'

nbDemos = size(Demos,2);
%% prepare GMM vs. GMM_delta figure
figure;
for ns = 1:4
    M = Gmms{ns};
    Mu_x = M.Mu;
    Sigma_x = M.Sigma;
    repro1 = M.repro;
    subplot(3,4,ns);hold on
    plotGMM(Mu_x([2,3],:), Sigma_x([2,3],[2,3],:), [0.5 .5 0.5], 1);
    for ii=1:nbDemos
        plot(Demos{ii}(1,:),Demos{ii}(2,:),'color',[0.5 0.5 0.5]);
    end
    
    plot(repro1(2,:),repro1(3,:),'r', 'linewidth',2);
    xticklabels([]);
    yticklabels([]);
    box on; grid on;
    title(['GMM #c' num2str(ns)],'fontname','Times','fontsize',14);
    subplot(3,4,[5,6,9,10]);hold on
    plot(repro1(2,:),repro1(3,:),'linewidth',2);
end

% return
subplot(3,4,[5,6,9,10]);hold on
for ii=1:nbDemos
    plot(Demos{ii}(1,:),Demos{ii}(2,:),'color',[0.5 0.5 0.5]);
end
xticklabels([]);
yticklabels([]);
box on; grid on;
title(['GMM '],'fontname','Times','fontsize',14);
subplot(3,4,[7,8,11,12]);hold on
for ii=1:size(Sols,2)
    plot(Sols{ii}(:,1),Sols{ii}(:,2),'linewidth',2);
end
for ii=1:nbDemos
    plot(Demos{ii}(1,:),Demos{ii}(2,:),'color',[0.5 0.5 0.5]);
end
xticklabels([]);
yticklabels([]);
box on; grid on;
title(['GMM-\delta '],'fontname','Times','fontsize',14);


subplot(3,4,1); hold on;
ylabel('x_2','fontname','Times','fontsize',14);
subplot(3,4,[5,6,9,10]); hold on;
ylabel('x_2','fontname','Times','fontsize',14);
xlabel('x_1','fontname','Times','fontsize',14);
subplot(3,4,[7,8,11,12]); hold on;
xlabel('x_1','fontname','Times','fontsize',14);

subplot(3,4,1);
axx = axis;
for ii=2:4
    subplot(3,4,ii);axis(axx);
end
subplot(3,4,[7,8,11,12]); axis auto; %axis(axx);
subplot(3,4,[5,6,9,10]);axis auto

% save fig file
return
fignamesaved = ['skill_' num2str(skillnum) '_compare_2Dw.fig'];
saveas(gcf,fignamesaved);

% print later for the paper