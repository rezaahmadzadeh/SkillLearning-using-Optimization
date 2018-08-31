function [Mu_g, R_Sigma_g_, G] = trainGMMG(Demos, nbDims, nbDemos, nbNodes, nbStates)
% Function to encode GMMs in gradient space of the demonstrations

G = -diag(ones(nbNodes,1)) + diag(ones(nbNodes-1,1),1);
% G = [G [zeros(nbNodes-2,1);1]];


g = cell(1,nbDemos);
for ii=1:nbDemos
    for jj=1:nbDims
        g{ii}(:,jj) = G * Demos{ii}(jj,:).';
    end
end

% Apply constraints (boundary condition for first/last sampling point)
% G(1,:)       = [];
G(end,:)     = [];
for ii=1:nbDemos
%     g{ii}(1,:) = [];
    g{ii}(end,:) = [];
end

% restructure data
% nbNodesMod = nbNodes - 1;
nbNodesMod = size(G,1);
D2 = zeros(nbDims+1, nbDemos*nbNodesMod);
% tt = linspace(0,1,nbNodesMod);
tt = 1:nbNodesMod;
D2(1,:) = repmat(tt,1,nbDemos);

for ii=1:nbDemos
    D2(2:nbDims+1,(ii-1)*nbNodesMod+1:ii*nbNodesMod) = g{ii}.'; %demos{ii}.pos;
end

% train GMM
M2 = encodeGMM(D2, nbNodesMod, nbStates);

% reproduce GMR
[repro2, expSigma2] = reproGMM(M2);

Mu_g = repro2(2:nbDims+1,:);
Sigma_g = expSigma2;
Sigma_g_ = [];
for ii=1:size(Sigma_g,3)
    Sigma_g_ = blkdiag(Sigma_g_,Sigma_g(:,:,ii));
end

R_Sigma_g_ = chol(inv(Sigma_g_));