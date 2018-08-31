function [Mu_d, R_Sigma_d_, L] = trainGMML(Demos, nbDims, nbDemos, nbNodes, nbStates)
L   = diag(ones(nbNodes,1)*2) - ...
    diag(ones(nbNodes-1,1),1) - ...
    diag(ones(nbNodes-1,1),-1);
L(1,2)       = -2;
L(end,end-1) = -2;
L            = L/2;

delta = cell(1,nbDemos);
for ii=1:nbDemos
    for jj=1:nbDims
        delta{ii}(:,jj) = L * Demos{ii}(jj,:).';
    end
end

% Apply constraints (boundary condition for first/last sampling point)
L(1,:)       = [];
L(end,:)     = [];
for ii=1:nbDemos
    delta{ii}(1,:) = [];
    delta{ii}(end,:) = [];
end

% restructure data
% nbNodesMod = nbNodes - 2;
nbNodesMod = size(L,1);
D2 = zeros(nbDims+1, nbDemos*nbNodesMod);
% tt = linspace(0,1,nbNodesMod);
tt = 1:nbNodesMod;
D2(1,:) = repmat(tt,1,nbDemos);

for ii=1:nbDemos
    D2(2:nbDims+1,(ii-1)*nbNodesMod+1:ii*nbNodesMod) = delta{ii}.'; %demos{ii}.pos;
end

% train GMM
M2 = encodeGMM(D2, nbNodesMod, nbStates);

% reproduce GMR
[repro2, expSigma2] = reproGMM(M2);

Mu_d = repro2(2:nbDims+1,:);
Sigma_d = expSigma2;
Sigma_d_ = [];
for ii=1:size(Sigma_d,3)
    Sigma_d_ = blkdiag(Sigma_d_,Sigma_d(:,:,ii));
end

R_Sigma_d_ = chol(inv(Sigma_d_));