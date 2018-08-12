function plotGMMGMR(data, numVar, M, expData, expSigma)
Mu = M.Mu;
Priors = M.Priors;
Sigma = M.Sigma;

%% Plot of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on
for n=1:numVar-1
  subplot(3*(numVar-1),2,(n-1)*2+1); hold on;
  plot(data(1,:), data(n+1,:), 'x', 'markerSize', 4, 'color', [.3 .3 .3]);
  axis([min(data(1,:)) max(data(1,:)) min(data(n+1,:))-0.01 max(data(n+1,:))+0.01]);
  xlabel('t','fontsize',16); ylabel(['x_' num2str(n)],'fontsize',16);
end
%plot 2D
subplot(3*(numVar-1),2,[2:2:2*(numVar-1)]); hold on;
plot(data(2,:), data(3,:), 'x', 'markerSize', 4, 'color', [.3 .3 .3]);
axis([min(data(2,:))-0.01 max(data(2,:))+0.01 min(data(3,:))-0.01 max(data(3,:))+0.01]);
xlabel('x_1','fontsize',16); ylabel('x_2','fontsize',16);

%% Plot of the GMM encoding results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot 1D
for n=1:numVar-1
  subplot(3*(numVar-1),2,4+(n-1)*2+1); hold on;
  plotGMM(Mu([1,n+1],:), Sigma([1,n+1],[1,n+1],:), [0 .8 0], 1);
  axis([min(data(1,:)) max(data(1,:)) min(data(n+1,:))-0.01 max(data(n+1,:))+0.01]);
  xlabel('t','fontsize',16); ylabel(['x_' num2str(n)],'fontsize',16);
end
%plot 2D
subplot(3*(numVar-1),2,4+[2:2:2*(numVar-1)]); hold on;
plotGMM(Mu([2,3],:), Sigma([2,3],[2,3],:), [0 .8 0], 1);
axis([min(data(2,:))-0.01 max(data(2,:))+0.01 min(data(3,:))-0.01 max(data(3,:))+0.01]);
xlabel('x_1','fontsize',16); ylabel('x_2','fontsize',16);

%% Plot of the GMR regression results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot 1D
for n=1:numVar-1
  subplot(3*(numVar-1),2,8+(n-1)*2+1); hold on;
  plotGMM(expData([1,n+1],:), expSigma(n,n,:), [0 0 .8], 3);
  axis([min(data(1,:)) max(data(1,:)) min(data(n+1,:))-0.01 max(data(n+1,:))+0.01]);
  xlabel('t','fontsize',16); ylabel(['x_' num2str(n)],'fontsize',16);
end
%plot 2D
subplot(3*(numVar-1),2,8+[2:2:2*(numVar-1)]); hold on;
plotGMM(expData([2,3],:), expSigma([1,2],[1,2],:), [0 0 .8], 2);
axis([min(data(2,:))-0.01 max(data(2,:))+0.01 min(data(3,:))-0.01 max(data(3,:))+0.01]);
xlabel('x_1','fontsize',16); ylabel('x_2','fontsize',16);

end