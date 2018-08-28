function [Mu_x, R_Sigma_x_] = trainGMM(Demos, nbDims, nbDemos, nbNodes, nbStates)
% restructure the data
D1 = zeros(nbDims+1, nbDemos*nbNodes);
t = 1:nbNodes;
D1(1,:) = repmat(t, 1, nbDemos);
for ii=1:nbDemos
    D1(2:nbDims+1, (ii-1)*nbNodes+1:ii*nbNodes) = Demos{ii};
end
% train GMM/GMR
M = encodeGMM(D1, nbNodes, nbStates);
[repro1, expSigma1] = reproGMM(M);
% Mu_x = M.Mu;
% % Priors = M.Priors;
% Sigma_x = M.Sigma;
%
% figure;hold on;
% plotGMM(Mu_x([2,3],:), Sigma_x([2,3],[2,3],:), [0 .8 0], 1);
% for ii=1:nbDemos
%     plot(Demos{ii}(1,:),Demos{ii}(2,:),'color',[0.5 0.5 0.5]);
% end
% plot(repro1(2,:),repro1(3,:),'r', 'linewidth',2);
Mu_x = repro1(2:nbDims+1,:);
Sigma_x = expSigma1;
Sigma_x_ = [];
for ii=1:size(Sigma_x,3)
    Sigma_x_ = blkdiag(Sigma_x_,Sigma_x(:,:,ii));
end

R_Sigma_x_ = chol(inv(Sigma_x_));
end